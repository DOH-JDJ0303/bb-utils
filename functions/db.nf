/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DATABASE SCANNING
    Reads the layout of a BigBacter database from file listings only:

        <db>/<taxa>/config.json
        <db>/<taxa>/sig/<sample>.sig.gz
        <db>/<taxa>/clusters/<cluster>/asm/<sample>.fa.gz
        <db>/<taxa>/clusters/<cluster>/asm/ref.fa.gz
        <db>/<taxa>/clusters/<cluster>/aux/ref.json.gz
        <db>/<taxa>/clusters/<cluster>/var/<sample>.tar.gz
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { listDir; runParallel; readJson } from './utils'

// Names of the taxa (top-level directories) in a database.
def listTaxa(db) {
    return listDir(db).findAll { info -> info.dir }.collect { info -> info.name }.sort()
}

// Strip a known suffix from a file name to get the sample name. Returns null if it does not match.
def sampleFromName(name, suffix) {
    if (!name.endsWith(suffix) || name.length() == suffix.length()) {
        return null
    }
    return name.substring(0, name.length() - suffix.length())
}

// Read the files of one cluster directory.
def scanCluster(clusterDir, clusterName) {
    def asm = listDir(clusterDir.resolve('asm')).findAll { info -> !info.dir }
    def var = listDir(clusterDir.resolve('var')).findAll { info -> !info.dir }
    def aux = listDir(clusterDir.resolve('aux')).findAll { info -> !info.dir }

    def asmBySample = [:]
    def varBySample = [:]
    def other = []

    asm.each { info ->
        def sample = sampleFromName(info.name, '.fa.gz')
        if (info.name == 'ref.fa.gz') {
            return
        }
        if (sample) {
            asmBySample[sample] = info
        }
        else {
            other << info
        }
    }
    var.each { info ->
        def sample = sampleFromName(info.name, '.tar.gz')
        if (sample) {
            varBySample[sample] = info
        }
        else {
            other << info
        }
    }

    // BigBacter writes ref.json.gz but also reads ref.json
    def refMeta = aux.find { info -> info.name == 'ref.json.gz' } ?: aux.find { info -> info.name == 'ref.json' }
    aux.each { info ->
        if (!(info.name in ['ref.json.gz', 'ref.json'])) {
            other << info
        }
    }

    return [
        name   : clusterName,
        path   : clusterDir,
        ref    : asm.find { info -> info.name == 'ref.fa.gz' },
        refMeta: refMeta,
        aux    : aux.findAll { info -> info.name in ['ref.json.gz', 'ref.json'] },
        asm    : asmBySample,
        var    : varBySample,
        other  : other,
        samples: (asmBySample.keySet() + varBySample.keySet()).toList().sort(),
    ]
}

// Read the reference metadata (name, length, n_contigs) of a scanned cluster, or null.
def readRefMeta(cluster) {
    if (!cluster.refMeta) {
        return null
    }
    try {
        return readJson(cluster.refMeta.path)
    }
    catch (Exception e) {
        return [error: e.message]
    }
}

// Scan one taxon. Returns its config, signatures, clusters and a flat list of sample records.
def scanTaxon(db, taxa, threads) {
    def dir = db.resolve(taxa)
    def top = listDir(dir)

    def configFile = top.find { info -> info.name == 'config.json' && !info.dir }
    def config = null
    def configError = null
    if (configFile) {
        try {
            config = readJson(configFile.path)
        }
        catch (Exception e) {
            configError = e.message
        }
    }

    def sigs = [:]
    listDir(dir.resolve('sig')).findAll { info -> !info.dir }.each { info ->
        def sample = sampleFromName(info.name, '.sig.gz')
        if (sample) {
            sigs[sample] = info
        }
    }

    def clusterDirs = listDir(dir.resolve('clusters')).findAll { info -> info.dir }.sort { info -> info.name }
    def scanned = runParallel(clusterDirs, threads) { info -> scanCluster(info.path, info.name) }
    def clusters = [:]
    scanned.each { c -> clusters[c.name] = c }

    // One record per sample per cluster. A sample's "added" time is the last-modified time of its
    // variant tarball (then assembly). The signature is NOT used: Floc rewrites every signature in
    // the taxon on each BigBacter push, so signature times reflect the latest push, not the first.
    def samples = []
    def clustersBySample = [:]
    clusters.each { name, c ->
        c.samples.each { sample ->
            def var = c.var[sample]
            def asm = c.asm[sample]
            def sig = sigs[sample]
            samples << [
                taxa   : taxa,
                cluster: name,
                sample : sample,
                var    : var,
                asm    : asm,
                sig    : sig,
                added  : var?.mtime ?: asm?.mtime ?: sig?.mtime,
                size   : (var?.size ?: 0L) + (asm?.size ?: 0L) + (sig?.size ?: 0L),
            ]
            if (!clustersBySample.containsKey(sample)) {
                clustersBySample[sample] = []
            }
            clustersBySample[sample] << name
        }
    }

    // Signatures with no cluster files
    def orphanSigs = sigs.findAll { sample, _info -> !clustersBySample.containsKey(sample) }

    return [
        taxa            : taxa,
        path            : dir,
        exists          : !top.isEmpty(),
        configFile      : configFile,
        config          : config,
        configError     : configError,
        sigs            : sigs,
        clusters        : clusters,
        samples         : samples,
        clustersBySample: clustersBySample,
        orphanSigs      : orphanSigs,
        otherTop        : top.findAll { info -> !(info.name in ['config.json', 'sig', 'clusters']) },
    ]
}

// Check that a requested taxon exists, with a helpful message if it does not.
def requireTaxon(db, taxa) {
    def available = listTaxa(db)
    if (!(taxa in available)) {
        def close = available.findAll { t -> t.toLowerCase() == taxa.toLowerCase() }
        def hint = close ? " Did you mean '${close[0]}'? Taxa names are case-sensitive." : ''
        error("Taxon '${taxa}' was not found in ${db}.${hint}\nAvailable taxa: ${available ? available.join(', ') : '(none)'}")
    }
    return taxa
}
