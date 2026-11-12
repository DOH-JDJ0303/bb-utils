/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODE: remove
    Removes samples from one taxon, selected by name, by cluster, or by when they were added.

    For each removed sample, its signature, assembly and variant tarball are deleted together.
    The signature is deleted first so that an interrupted run never leaves a sample that Floc
    can still assign to but that has no cluster files.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { cleanTaxa; asList; parseTime; fmtTime; fmtSize; listRecursive; runParallel; writeTsv; textTable } from './utils'
include { scanTaxon; readRefMeta; requireTaxon                                                           } from './db'

// Sample names from --samples: a comma-separated list, or a file with one name per line.
def readSampleNames(value) {
    if (value == null || value.toString().trim() == '') {
        return null
    }
    def text = value.toString().trim()
    def names = null
    if (!text.contains(',')) {
        def candidate = file(text)
        if (candidate.exists() && candidate.isFile()) {
            names = candidate.readLines().collect { l -> l.trim() }.findAll { l -> l && !l.startsWith('#') }
        }
    }
    if (names == null) {
        names = asList(text)
    }
    if (!names) {
        error('--samples was given but no sample names were found.')
    }
    return names.unique()
}

def deleteOne(row) {
    try {
        java.nio.file.Files.deleteIfExists(row.pathObj)
        return 'deleted'
    }
    catch (Exception e) {
        return "failed: ${e.message}"
    }
}

def runRemove(opts) {
    def db = opts.db

    // ---------------------------------------------------------------------------------
    //  Validate the request
    // ---------------------------------------------------------------------------------
    if (!opts.taxa) {
        error('--mode remove needs --taxa. Removals are always limited to one taxon.')
    }
    if (opts.taxa.toString().contains(',')) {
        error('--mode remove works on one taxon at a time. Give a single --taxa.')
    }
    def taxa = requireTaxon(db, cleanTaxa(opts.taxa))
    def cluster = opts.cluster == null || opts.cluster.toString().trim() == '' ? null : opts.cluster.toString().trim()
    def names = readSampleNames(opts.samples)
    def after = parseTime(opts.added_after, 'added_after')
    def before = parseTime(opts.added_before, 'added_before')
    def byTime = after != null || before != null

    if (names && !cluster) {
        error('--samples needs --cluster. Samples are removed by name from one taxon and cluster.')
    }
    if (after != null && before != null && after >= before) {
        error("--added_after (${fmtTime(after)}) must be earlier than --added_before (${fmtTime(before)}).")
    }
    if (opts.all && (cluster || names || byTime)) {
        error('--all removes the whole taxon and cannot be combined with --cluster, --samples, --added_after or --added_before.')
    }
    if (!opts.all && !cluster && !names && !byTime) {
        error("Nothing selected. With --taxa ${taxa}, add one of:\n" +
              "  --cluster <name>                          remove a whole cluster\n" +
              "  --cluster <name> --samples <a,b,c|file>   remove samples by name\n" +
              "  --added_after / --added_before <time>     remove samples by when they were added\n" +
              "  --all                                     remove the whole taxon")
    }

    log.info("Scanning ${taxa} ...")
    def scan = scanTaxon(db, taxa, opts.threads)
    if (cluster && !scan.clusters.containsKey(cluster)) {
        error("Cluster '${cluster}' was not found in ${taxa}. Available clusters: ${scan.clusters.keySet().join(', ') ?: '(none)'}")
    }

    def rows = []
    def addRow = { clusterName, sample, type, info ->
        rows << [
            action       : 'delete',
            taxa         : taxa,
            cluster      : clusterName ?: '',
            sample       : sample ?: '',
            type         : type,
            path         : info.path.toUriString(),
            size_bytes   : info.size,
            last_modified: fmtTime(info.mtime),
            status       : 'planned',
            pathObj      : info.path,
        ]
    }
    def notes = []
    def selected = []
    def pruned = []
    def notFound = []

    if (opts.all) {
        // -----------------------------------------------------------------------------
        //  Whole taxon
        // -----------------------------------------------------------------------------
        selected = scan.samples
        pruned = scan.clusters.keySet().toList()
        listRecursive(scan.path).each { info -> addRow.call('', '', 'taxon', info) }
    }
    else {
        // -----------------------------------------------------------------------------
        //  Select samples
        // -----------------------------------------------------------------------------
        selected = scan.samples.findAll { s -> cluster == null || s.cluster == cluster }
        if (names) {
            notFound = names.findAll { n -> !selected.any { s -> s.sample == n } }
            selected = selected.findAll { s -> s.sample in names }
        }
        if (byTime) {
            selected = selected.findAll { s ->
                s.added != null && (after == null || s.added >= after) && (before == null || s.added < before)
            }
        }
        if (notFound) {
            notes << "WARNING: ${notFound.size()} sample(s) were not found in ${taxa} cluster ${cluster}: ${notFound.join(', ')}"
        }

        def selectedKeys = selected.collect { s -> "${s.cluster}\t${s.sample}".toString() }.toSet()
        def affected = selected.collect { s -> s.cluster }.unique().sort()

        // A signature is only deleted if the sample does not remain in another cluster
        def sigDone = [] as Set
        selected.each { s ->
            def stays = scan.clustersBySample[s.sample].any { c -> !selectedKeys.contains("${c}\t${s.sample}".toString()) }
            if (s.sig && !stays && !sigDone.contains(s.sample)) {
                addRow.call(s.cluster, s.sample, 'sig', s.sig)
                sigDone << s.sample
            }
            if (s.sig && stays) {
                notes << "NOTE: signature of ${s.sample} is kept because the sample is also in another cluster."
            }
        }
        selected.each { s ->
            if (s.var) {
                addRow.call(s.cluster, s.sample, 'var', s.var)
            }
            if (s.asm) {
                addRow.call(s.cluster, s.sample, 'asm', s.asm)
            }
        }

        // -----------------------------------------------------------------------------
        //  Clusters left with no samples
        // -----------------------------------------------------------------------------
        affected.each { name ->
            def c = scan.clusters[name]
            def remaining = c.samples.findAll { sample -> !selectedKeys.contains("${name}\t${sample}".toString()) }
            if (!remaining && opts.prune_empty) {
                pruned << name
                def already = rows.collect { r -> r.path }.toSet()
                listRecursive(c.path).each { info ->
                    if (!already.contains(info.path.toUriString())) {
                        def type = info.name == 'ref.fa.gz' ? 'ref' : (info.name in ['ref.json.gz', 'ref.json'] ? 'aux' : 'other')
                        addRow.call(name, '', type, info)
                    }
                }
            }
            else if (!remaining) {
                notes << "WARNING: cluster ${name} will have no samples left but is kept (--prune_empty false). Floc may reuse its name for unrelated samples, which would then inherit this cluster's reference."
            }
            else {
                // Removing the sample the reference was made from is allowed: ref.fa.gz is a separate, immutable copy
                def meta = readRefMeta(c)
                def refName = meta?.name
                if (refName && selectedKeys.contains("${name}\t${refName}".toString())) {
                    notes << "NOTE: ${refName} is the reference for cluster ${name}. Its sample files are removed, but ref.fa.gz is kept so the cluster stays usable."
                }
            }
        }
    }

    // ---------------------------------------------------------------------------------
    //  Report the plan and write the manifest BEFORE anything is deleted
    // ---------------------------------------------------------------------------------
    def columns = ['action', 'taxa', 'cluster', 'sample', 'type', 'path', 'size_bytes', 'last_modified', 'status']
    def manifest = opts.outdir.resolve('remove_manifest.tsv')
    def totalSize = rows.collect { r -> r.size_bytes }.sum() ?: 0L

    def lines = []
    lines << ''
    lines << "Database : ${db}"
    lines << "Taxon    : ${taxa}"
    if (opts.all) {
        lines << 'Selection: the whole taxon (--all)'
    }
    else {
        def how = []
        if (cluster) { how << "cluster ${cluster}" }
        if (names) { how << "${names.size()} named sample(s)" }
        if (after != null) { how << "added at or after ${fmtTime(after)}" }
        if (before != null) { how << "added before ${fmtTime(before)}" }
        lines << "Selection: ${how.join(', ')}"
    }
    lines << ''
    if (selected) {
        def perCluster = selected.groupBy { s -> s.cluster }.collect { name, list ->
            def added = list.collect { s -> s.added }.findAll { t -> t != null }
            [
                cluster      : name,
                remove       : list.size(),
                of           : scan.clusters[name].samples.size(),
                first_added  : fmtTime(added ? added.min() : null),
                last_added   : fmtTime(added ? added.max() : null),
                whole_cluster: name in pruned ? 'yes (cluster deleted)' : '',
            ]
        }.sort { r -> r.cluster.toString().padLeft(12, '0') }
        lines << textTable(['cluster', 'remove', 'of', 'first_added', 'last_added', 'whole_cluster'], perCluster)
        lines << ''
    }
    lines << "${selected.size()} sample(s) in ${selected.collect { s -> s.cluster }.unique().size()} cluster(s); ${pruned.size()} cluster(s) deleted entirely; ${rows.size()} file(s), ${fmtSize(totalSize)}"
    notes.unique().each { n -> lines << n }
    log.info(lines.join('\n'))

    if (!rows) {
        writeTsv(manifest, columns, rows)
        log.warn("No files matched the selection. Nothing was removed. (${manifest})")
        return
    }

    if (opts.dry_run) {
        rows.each { r -> r.status = 'dry_run' }
        writeTsv(manifest, columns, rows)
        log.info("\nDRY RUN: nothing was deleted. Files that would be deleted are listed in ${manifest}\n")
        return
    }
    writeTsv(manifest, columns, rows)

    // ---------------------------------------------------------------------------------
    //  Delete: signatures first, then sample files, then what is left of emptied clusters
    // ---------------------------------------------------------------------------------
    log.info("\nDeleting ${rows.size()} file(s) ...")
    def phases = [
        rows.findAll { r -> r.type == 'sig' },
        rows.findAll { r -> r.type in ['var', 'asm'] },
        rows.findAll { r -> !(r.type in ['sig', 'var', 'asm']) },
    ]
    phases.each { phase ->
        def results = runParallel(phase, opts.threads) { r -> deleteOne(r) }
        phase.eachWithIndex { r, i -> r.status = results[i] }
    }
    writeTsv(manifest, columns, rows)

    // Tidy up empty directories (only meaningful on file systems that have real directories)
    def failed = rows.findAll { r -> r.status != 'deleted' }
    if (!failed) {
        try {
            if (opts.all) {
                scan.path.deleteDir()
            }
            else {
                pruned.each { name -> scan.clusters[name].path.deleteDir() }
            }
        }
        catch (Exception e) {
            log.warn("Files were deleted, but empty directories could not be tidied up: ${e.message}")
        }
    }

    if (failed) {
        error("${failed.size()} of ${rows.size()} file(s) could not be deleted. See the status column in ${manifest}")
    }
    log.info("Deleted ${rows.size()} file(s), ${fmtSize(totalSize)}. Manifest: ${manifest}\n")
}
