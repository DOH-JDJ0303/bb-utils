/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MODE: summary
    Describes a BigBacter database and checks it for problems. Read-only.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { cleanTaxa; asList; fmtTime; fmtSize; runParallel; writeTsv; textTable } from './utils'
include { listTaxa; scanTaxon; readRefMeta; requireTaxon                         } from './db'

def runSummary(opts) {
    def db = opts.db
    def taxaList = listTaxa(db)
    if (opts.taxa) {
        taxaList = asList(opts.taxa).collect { t -> requireTaxon(db, cleanTaxa(t)) }
    }
    if (!taxaList) {
        error("No taxa found in ${db}. Check the --db path.")
    }
    def clusterFilter = opts.cluster == null ? null : asList(opts.cluster)

    def taxaRows = []
    def clusterRows = []
    def sampleRows = []
    def issues = []

    taxaList.each { taxa ->
        log.info("Scanning ${taxa} ...")
        def scan = scanTaxon(db, taxa, opts.threads)
        def clusters = scan.clusters.findAll { name, _c -> clusterFilter == null || name in clusterFilter }
        def samples = scan.samples.findAll { s -> clusters.containsKey(s.cluster) }

        // Reference metadata: one small read per cluster
        def refMeta = [:]
        def clusterList = clusters.values().toList()
        def metas = runParallel(clusterList, opts.threads) { c -> readRefMeta(c) }
        clusterList.eachWithIndex { c, i -> refMeta[c.name] = metas[i] }

        // ---- taxon-level checks
        def taxonIssue = { msg -> issues << [level: 'taxon', taxa: taxa, cluster: '', sample: '', issue: msg] }
        if (!scan.configFile && (scan.sigs || scan.clusters)) {
            taxonIssue.call('config.json is missing: BigBacter v2.2.0+ will not run against this taxon')
        }
        else if (scan.configError) {
            taxonIssue.call("config.json could not be read: ${scan.configError}")
        }
        else if (scan.configFile && scan.config?.clust_dist == null) {
            taxonIssue.call('config.json does not record clust_dist')
        }
        if (clusterFilter == null) {
            scan.orphanSigs.each { sample, _info ->
                issues << [level: 'sample', taxa: taxa, cluster: '', sample: sample, issue: 'signature has no assembly or variant files in any cluster']
            }
        }

        // ---- cluster-level rows and checks
        clusters.each { name, c ->
            def members = samples.findAll { s -> s.cluster == name }
            def meta = refMeta[name]
            def clusterIssue = { msg -> issues << [level: 'cluster', taxa: taxa, cluster: name, sample: '', issue: msg] }
            if (!c.ref) {
                clusterIssue.call('ref.fa.gz is missing: BigBacter will fail for this cluster')
            }
            if (!c.refMeta) {
                clusterIssue.call('ref.json.gz is missing')
            }
            else if (meta?.error) {
                clusterIssue.call("ref.json could not be read: ${meta.error}")
            }
            if (!members) {
                clusterIssue.call('cluster has no samples: its name could be reused for unrelated samples')
            }
            def added = members.collect { s -> s.added }.findAll { t -> t != null }
            def refSize = (c.ref?.size ?: 0L) + (c.aux.collect { a -> a.size }.sum() ?: 0L) + (c.other.collect { a -> a.size }.sum() ?: 0L)
            def size = refSize + (members.collect { s -> (s.var?.size ?: 0L) + (s.asm?.size ?: 0L) }.sum() ?: 0L)
            clusterRows << [
                taxa        : taxa,
                cluster     : name,
                samples     : members.size(),
                ref_name    : meta?.name,
                ref_length  : meta?.length,
                ref_contigs : meta?.n_contigs,
                size_bytes  : size,
                size        : fmtSize(size),
                first_added : fmtTime(added ? added.min() : null),
                last_added  : fmtTime(added ? added.max() : null),
            ]
        }

        // ---- sample-level rows and checks
        samples.each { s ->
            def sampleIssue = { msg -> issues << [level: 'sample', taxa: taxa, cluster: s.cluster, sample: s.sample, issue: msg] }
            if (!s.var) {
                sampleIssue.call('variant tarball (var) is missing')
            }
            if (!s.asm) {
                sampleIssue.call('assembly (asm) is missing')
            }
            if (!s.sig) {
                sampleIssue.call('signature (sig) is missing: expected only when the cluster was set in the samplesheet')
            }
            if (scan.clustersBySample[s.sample].size() > 1) {
                sampleIssue.call("sample is in more than one cluster: ${scan.clustersBySample[s.sample].join(', ')}")
            }
            sampleRows << [
                taxa      : taxa,
                cluster   : s.cluster,
                sample    : s.sample,
                added     : fmtTime(s.added),
                has_sig   : s.sig != null,
                has_asm   : s.asm != null,
                has_var   : s.var != null,
                size_bytes: s.size,
            ]
        }

        // ---- taxon-level row
        def added = samples.collect { s -> s.added }.findAll { t -> t != null }
        def sigSize = clusterFilter == null ? (scan.sigs.values().collect { i -> i.size }.sum() ?: 0L) : (samples.collect { s -> s.sig?.size ?: 0L }.sum() ?: 0L)
        def size = (clusterRows.findAll { r -> r.taxa == taxa }.collect { r -> r.size_bytes }.sum() ?: 0L) + sigSize + (scan.configFile?.size ?: 0L)
        taxaRows << [
            taxa                 : taxa,
            clusters             : clusters.size(),
            samples              : samples.collect { s -> s.sample }.unique().size(),
            signatures           : clusterFilter == null ? scan.sigs.size() : samples.count { s -> s.sig != null },
            size_bytes           : size,
            size                 : fmtSize(size),
            first_added          : fmtTime(added ? added.min() : null),
            last_added           : fmtTime(added ? added.max() : null),
            clust_dist           : scan.config?.clust_dist,
            clust_assign_hashfreq: scan.config?.clust_assign_hashfreq,
            floc_version         : scan.config?.floc_version,
            bigbacter_version    : scan.config?.bigbacter_version,
            issues               : issues.count { i -> i.taxa == taxa },
        ]
    }

    // ---- write tables
    def out = opts.outdir
    def taxaCols    = ['taxa', 'clusters', 'samples', 'signatures', 'size_bytes', 'size', 'first_added', 'last_added', 'clust_dist', 'clust_assign_hashfreq', 'floc_version', 'bigbacter_version', 'issues']
    def clusterCols = ['taxa', 'cluster', 'samples', 'ref_name', 'ref_length', 'ref_contigs', 'size_bytes', 'size', 'first_added', 'last_added']
    def sampleCols  = ['taxa', 'cluster', 'sample', 'added', 'has_sig', 'has_asm', 'has_var', 'size_bytes']
    def issueCols   = ['level', 'taxa', 'cluster', 'sample', 'issue']
    writeTsv(out.resolve('summary_taxa.tsv'), taxaCols, taxaRows)
    writeTsv(out.resolve('summary_clusters.tsv'), clusterCols, clusterRows)
    writeTsv(out.resolve('summary_samples.tsv'), sampleCols, sampleRows)
    writeTsv(out.resolve('summary_issues.tsv'), issueCols, issues)

    // ---- console report
    def total = taxaRows.collect { r -> r.size_bytes }.sum() ?: 0L
    def lines = []
    lines << ''
    lines << "BigBacter database: ${db}"
    lines << "  ${taxaRows.size()} taxa, ${clusterRows.size()} clusters, ${sampleRows.size()} samples, ${fmtSize(total)}"
    lines << ''
    lines << textTable(['taxa', 'clusters', 'samples', 'size', 'first_added', 'last_added', 'clust_dist', 'issues'], taxaRows)
    lines << ''
    if (issues) {
        lines << "${issues.size()} issue(s) found:"
        def grouped = issues.groupBy { i -> i.issue.replaceAll(':.*', '') }
        grouped.each { msg, list -> lines << "  ${list.size().toString().padLeft(6)}  ${msg}" }
        lines << "  See ${out.resolve('summary_issues.tsv')}"
    }
    else {
        lines << 'No issues found.'
    }
    lines << ''
    lines << "Tables written to ${out}"
    log.info(lines.join('\n'))
}
