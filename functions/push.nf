/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODE: push
    Copies the database files of one BigBacter run (<outdir>/<timestamp>) into the database.

    Run directory                                   Database
    ----------------------------------------------  -------------------------------------------------
    <taxa>/config.json                              <taxa>/config.json
    <taxa>/sig/<sample>.sig.gz                      <taxa>/sig/<sample>.sig.gz
    <taxa>/<cluster>/asm/ref.fa.gz                  <taxa>/clusters/<cluster>/asm/ref.fa.gz
    <taxa>/<cluster>/aux/ref.json.gz                <taxa>/clusters/<cluster>/aux/ref.json.gz
    <taxa>/<cluster>/var/<sample>.tar.gz            <taxa>/clusters/<cluster>/var/<sample>.tar.gz
    <taxa>/<cluster>/asm/<sample>.fa.gz             <taxa>/clusters/<cluster>/asm/<sample>.fa.gz

    The samples to push are the variant tarballs found in the run: BigBacter only publishes a
    tarball for samples analysed in that run, so no summary file is needed to find them.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { cleanTaxa; asList; fmtSize; listDir; runParallel; readJson; sameNumber; writeTsv; textTable } from './utils'
include { scanTaxon; sampleFromName                                                                                  } from './db'

// Directories directly inside a run directory that are not taxa
def sampleDirNames() {
    return ['asm', 'taxa', 'reads'] as Set
}

// Read the database files that one run holds for one taxon.
def scanRunTaxon(taxaInfo, threads) {
    def entries = listDir(taxaInfo.path)
    def sigs = [:]
    listDir(taxaInfo.path.resolve('sig')).findAll { info -> !info.dir }.each { info ->
        def sample = sampleFromName(info.name, '.sig.gz')
        if (sample) {
            sigs[sample] = info
        }
    }

    def candidates = entries.findAll { info -> info.dir && !(info.name in ['sig', 'cluster']) }
    def scanned = runParallel(candidates, threads) { info ->
        def var = [:]
        listDir(info.path.resolve('var')).findAll { f -> !f.dir }.each { f ->
            def sample = sampleFromName(f.name, '.tar.gz')
            if (sample && sample != 'ref') {
                var[sample] = f
            }
        }
        def asm = listDir(info.path.resolve('asm'))
        def aux = listDir(info.path.resolve('aux'))
        // Sample assemblies sit next to the reference: <taxa>/<cluster>/asm/<sample>.fa.gz
        def sampleAsm = [:]
        asm.findAll { f -> !f.dir && f.name != 'ref.fa.gz' }.each { f ->
            def sample = sampleFromName(f.name, '.fa.gz')
            if (sample) {
                sampleAsm[sample] = f
            }
        }
        return [
            name   : info.name,
            var    : var,
            asm    : sampleAsm,
            ref    : asm.find { f -> f.name == 'ref.fa.gz' },
            refMeta: aux.find { f -> f.name == 'ref.json.gz' } ?: aux.find { f -> f.name == 'ref.json' },
        ]
    }
    def clusters = scanned.findAll { c -> c.var }

    return [
        taxa      : taxaInfo.name,
        configFile: entries.find { info -> info.name == 'config.json' && !info.dir },
        sigs      : sigs,
        clusters  : clusters,
    ]
}

def copyOne(row) {
    try {
        def target = row.destObj
        if (target.toUri().getScheme() == 'file') {
            java.nio.file.Files.createDirectories(target.getParent())
        }
        row.sourceObj.copyTo(target)
        return 'copied'
    }
    catch (Exception e) {
        return "failed: ${e.message}"
    }
}

def runPush(opts) {
    def db = opts.db
    if (!opts.run) {
        error('--mode push needs --run: the BigBacter run directory to push, i.e. <outdir>/<timestamp>.')
    }
    def runDir = opts.run
    if (!runDir.exists() || !runDir.isDirectory()) {
        error("--run ${runDir.toUriString()} is not a directory.")
    }

    // ---------------------------------------------------------------------------------
    //  Find the taxa in the run
    // ---------------------------------------------------------------------------------
    def top = listDir(runDir).findAll { info -> info.dir }
    def wanted = opts.taxa ? asList(opts.taxa).collect { t -> cleanTaxa(t) } : null
    if (wanted) {
        def missing = wanted.findAll { t -> !top.any { info -> info.name == t } }
        if (missing) {
            error("--taxa ${missing.join(', ')} not found in ${runDir.toUriString()}.")
        }
        top = top.findAll { info -> info.name in wanted }
    }
    else {
        // Skip the per-sample folders (they only hold taxa/ and reads/, and asm/ in runs before this layout)
        def kinds = runParallel(top, opts.threads) { info ->
            def children = listDir(info.path).collect { c -> c.name }
            return children && !sampleDirNames().containsAll(children)
        }
        def keep = []
        top.eachWithIndex { info, i -> if (kinds[i]) { keep << info } }
        top = keep
    }

    def rows = []
    def addRow = { taxa, cluster, sample, type, source, dest, status ->
        rows << [
            action    : status == 'planned' ? 'copy' : 'skip',
            taxa      : taxa,
            cluster   : cluster ?: '',
            sample    : sample ?: '',
            type      : type,
            source    : source ? source.path.toUriString() : '',
            dest      : dest ? dest.toUriString() : '',
            size_bytes: source ? source.size : '',
            status    : status,
            sourceObj : source?.path,
            destObj   : dest,
        ]
    }
    def problems = []
    def warnings = []
    def overview = []

    top.each { taxaInfo ->
        def taxa = taxaInfo.name
        def run = scanRunTaxon(taxaInfo, opts.threads)
        if (!run.clusters) {
            return
        }
        log.info("Checking ${taxa} ...")
        def dbScan = scanTaxon(db, taxa, opts.threads)
        def dbTaxon = db.resolve(taxa)
        def nSamples = run.clusters.collect { c -> c.var.size() }.sum()
        def stats = [taxa: taxa, clusters: run.clusters.size(), samples: nSamples, push: 0, in_db: 0, skipped: 0, new_clusters: 0]
        overview << stats

        // -----------------------------------------------------------------------------
        //  config.json: same rules as BigBacter
        // -----------------------------------------------------------------------------
        def runConfig = null
        if (run.configFile) {
            try {
                runConfig = readJson(run.configFile.path)
            }
            catch (Exception e) {
                problems << "${taxa}: config.json in the run could not be read (${e.message}). Taxon skipped."
                stats.skipped = nSamples
                return
            }
        }
        def dbHasData = dbScan.sigs || dbScan.clusters
        def taxonError = null
        if (dbScan.configFile) {
            if (dbScan.configError || dbScan.config?.clust_dist == null) {
                taxonError = 'config.json in the database could not be read or has no clust_dist'
            }
            else if (runConfig && !sameNumber(runConfig.clust_dist, dbScan.config.clust_dist)) {
                taxonError = "the run used clust_dist ${runConfig.clust_dist} but the database was built with ${dbScan.config.clust_dist}"
            }
            else if (runConfig && runConfig.clust_assign_hashfreq != null && dbScan.config.clust_assign_hashfreq != null && !sameNumber(runConfig.clust_assign_hashfreq, dbScan.config.clust_assign_hashfreq)) {
                warnings << "${taxa}: the run used clust_assign_hashfreq ${runConfig.clust_assign_hashfreq} but the database was built with ${dbScan.config.clust_assign_hashfreq}."
            }
        }
        else if (dbHasData) {
            taxonError = 'the database already has this taxon but it has no config.json (see the BigBacter v2.2.0 upgrade notes)'
        }
        else if (!run.configFile) {
            taxonError = 'the taxon is new to the database and the run has no config.json'
        }
        if (taxonError) {
            problems << "${taxa}: ${taxonError}. Taxon skipped."
            stats.skipped = nSamples
            run.clusters.each { c -> c.var.each { sample, info -> addRow.call(taxa, c.name, sample, 'var', info, null, "skipped: ${taxonError}") } }
            return
        }
        if (!dbScan.configFile) {
            addRow.call(taxa, '', '', 'config', run.configFile, dbTaxon.resolve('config.json'), 'planned')
        }

        run.clusters.each { c ->
            def dbCluster = dbScan.clusters[c.name]
            def dest = dbTaxon.resolve('clusters').resolve(c.name)

            // -------------------------------------------------------------------------
            //  Reference: never overwritten, and must be the same one the run used
            // -------------------------------------------------------------------------
            def clusterError = null
            if (dbCluster?.ref) {
                if (c.ref && c.ref.size != dbCluster.ref.size) {
                    clusterError = 'the reference in the run differs from the one in the database (was the cluster name reused?)'
                }
            }
            else if (dbCluster && dbCluster.samples) {
                clusterError = 'the cluster exists in the database with samples but has no ref.fa.gz'
            }
            else if (!c.ref) {
                clusterError = 'the cluster is new to the database and the run has no ref.fa.gz for it'
            }
            if (clusterError) {
                problems << "${taxa} cluster ${c.name}: ${clusterError}. Cluster skipped (${c.var.size()} sample(s))."
                stats.skipped = stats.skipped + c.var.size()
                c.var.each { sample, info -> addRow.call(taxa, c.name, sample, 'var', info, null, "skipped: ${clusterError}") }
                return
            }
            if (!dbCluster?.ref) {
                stats.new_clusters = stats.new_clusters + 1
                addRow.call(taxa, c.name, '', 'ref', c.ref, dest.resolve('asm').resolve('ref.fa.gz'), 'planned')
            }
            if (!dbCluster?.refMeta) {
                if (c.refMeta) {
                    addRow.call(taxa, c.name, '', 'aux', c.refMeta, dest.resolve('aux').resolve(c.refMeta.name), 'planned')
                }
                else {
                    warnings << "${taxa} cluster ${c.name}: no ref.json.gz in the run or the database. Reports will label the reference row 'Reference'."
                }
            }

            // -------------------------------------------------------------------------
            //  Samples
            // -------------------------------------------------------------------------
            c.var.each { sample, varInfo ->
                def asmInfo = c.asm[sample]
                def sigInfo = run.sigs[sample]
                def elsewhere = (dbScan.clustersBySample[sample] ?: []).findAll { name -> name != c.name }
                def sampleError = null
                if (elsewhere) {
                    sampleError = "already in the database under cluster ${elsewhere.join(', ')}"
                }
                else if (!asmInfo) {
                    sampleError = "assembly not found in the run (${taxa}/${c.name}/asm/${sample}.fa.gz)"
                }
                else if (!sigInfo && !opts.allow_missing_sig) {
                    sampleError = "signature not found in the run (${taxa}/sig/${sample}.sig.gz). Use --allow_missing_sig for samples whose cluster was set in the samplesheet"
                }
                if (sampleError) {
                    problems << "${taxa} cluster ${c.name} sample ${sample}: ${sampleError}. Sample skipped."
                    stats.skipped = stats.skipped + 1
                    addRow.call(taxa, c.name, sample, 'var', varInfo, null, "skipped: ${sampleError}")
                    return
                }
                if (!sigInfo) {
                    warnings << "${taxa} cluster ${c.name} sample ${sample}: pushed without a signature."
                }

                def parts = [
                    [type: 'var', source: varInfo, dest: dest.resolve('var').resolve(varInfo.name), inDb: dbCluster?.var?.get(sample) != null],
                    [type: 'asm', source: asmInfo, dest: dest.resolve('asm').resolve(asmInfo.name), inDb: dbCluster?.asm?.get(sample) != null],
                ]
                if (sigInfo) {
                    parts << [type: 'sig', source: sigInfo, dest: dbTaxon.resolve('sig').resolve(sigInfo.name), inDb: dbScan.sigs[sample] != null]
                }
                def todo = parts.findAll { p -> !p.inDb || opts.overwrite }
                parts.each { p ->
                    def status = (!p.inDb || opts.overwrite) ? 'planned' : 'skipped: already in the database'
                    addRow.call(taxa, c.name, sample, p.type, p.source, p.dest, status)
                }
                if (todo) {
                    stats.push = stats.push + 1
                }
                else {
                    stats.in_db = stats.in_db + 1
                }
            }
        }
    }

    // ---------------------------------------------------------------------------------
    //  Report the plan and write the manifest BEFORE anything is copied
    // ---------------------------------------------------------------------------------
    def columns = ['action', 'taxa', 'cluster', 'sample', 'type', 'source', 'dest', 'size_bytes', 'status']
    def manifest = opts.outdir.resolve('push_manifest.tsv')
    def todo = rows.findAll { r -> r.status == 'planned' }
    def totalSize = todo.collect { r -> r.size_bytes }.sum() ?: 0L

    def lines = []
    lines << ''
    lines << "Run      : ${runDir.toUriString()}"
    lines << "Database : ${db.toUriString()}"
    lines << ''
    if (overview) {
        lines << textTable(['taxa', 'clusters', 'new_clusters', 'samples', 'push', 'in_db', 'skipped'], overview)
        lines << ''
    }
    lines << "${todo.size()} file(s) to copy, ${fmtSize(totalSize)}" + (opts.overwrite ? ' (--overwrite: existing sample files are replaced)' : '')
    warnings.each { w -> lines << "WARNING: ${w}" }
    problems.each { p -> lines << "PROBLEM: ${p}" }
    log.info(lines.join('\n'))

    if (!overview) {
        error("No database files were found in ${runDir.toUriString()}. --run should point at one BigBacter run, i.e. <outdir>/<timestamp>, holding <taxa>/<cluster>/var/<sample>.tar.gz.")
    }

    if (opts.dry_run) {
        todo.each { r -> r.status = 'dry_run' }
        writeTsv(manifest, columns, rows)
        log.info("\nDRY RUN: nothing was copied. The plan is listed in ${manifest}\n")
        return
    }
    writeTsv(manifest, columns, rows)

    // ---------------------------------------------------------------------------------
    //  Copy: config and references first, then sample files, signatures last.
    //  A sample only becomes visible to Floc (through its signature) once its other files exist.
    // ---------------------------------------------------------------------------------
    if (todo) {
        log.info("\nCopying ${todo.size()} file(s) ...")
    }
    def runPhase = { phase ->
        def results = runParallel(phase, opts.threads) { r -> copyOne(r) }
        phase.eachWithIndex { r, i -> r.status = results[i] }
    }
    def key = { r -> "${r.taxa}\t${r.cluster}".toString() }

    def configRows = todo.findAll { r -> r.type == 'config' }
    runPhase.call(configRows)
    def badTaxa = configRows.findAll { r -> r.status != 'copied' }.collect { r -> r.taxa }.toSet()

    def refRows = todo.findAll { r -> r.type in ['ref', 'aux'] && !badTaxa.contains(r.taxa) }
    runPhase.call(refRows)
    def badClusters = refRows.findAll { r -> r.type == 'ref' && r.status != 'copied' }.collect { r -> key.call(r) }.toSet()

    def blocked = { r -> badTaxa.contains(r.taxa) || badClusters.contains(key.call(r)) }
    def fileRows = todo.findAll { r -> r.type in ['var', 'asm'] && !blocked.call(r) }
    runPhase.call(fileRows)
    def badSamples = fileRows.findAll { r -> r.status != 'copied' }.collect { r -> "${key.call(r)}\t${r.sample}".toString() }.toSet()

    def sigRows = todo.findAll { r -> r.type == 'sig' && !blocked.call(r) && !badSamples.contains("${key.call(r)}\t${r.sample}".toString()) }
    runPhase.call(sigRows)

    todo.findAll { r -> r.status == 'planned' }.each { r -> r.status = 'not copied: an earlier file for this taxon, cluster or sample failed' }
    writeTsv(manifest, columns, rows)

    def copied = todo.findAll { r -> r.status == 'copied' }
    def failed = todo.findAll { r -> r.status != 'copied' }
    log.info("Copied ${copied.size()} file(s), ${fmtSize(copied.collect { r -> r.size_bytes }.sum() ?: 0L)}. Manifest: ${manifest}\n")
    if (failed) {
        error("${failed.size()} of ${todo.size()} file(s) could not be copied. See the status column in ${manifest}")
    }
    if (problems) {
        error("The push finished, but ${problems.size()} taxon, cluster or sample problem(s) were skipped (listed above and in ${manifest}).")
    }
}
