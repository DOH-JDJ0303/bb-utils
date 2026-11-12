/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENERAL HELPERS
    Everything here uses the Java NIO file API through Nextflow, so the same code works
    on local paths, s3://, gs:// and az://. Nothing is staged or downloaded: listing and
    deleting only touch object metadata.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Clean a taxa name the same way BigBacter does for samplesheet and GAMBIT names
// (subworkflows/local/utils_nfcore_bigbacter_pipeline and subworkflows/local/classify).
def cleanTaxa(taxa) {
    if (taxa == null) {
        return null
    }
    return taxa.toString().trim().replaceAll('[^a-zA-Z0-9]', '_')
}

// Interpret a parameter as true/false, whether it arrives as a Boolean or a String.
def asBool(value) {
    if (value == null) {
        return false
    }
    if (value instanceof Boolean) {
        return value
    }
    return value.toString().trim().toLowerCase() in ['true', 'yes', '1']
}

// Split a comma-separated parameter into a list of trimmed, non-empty strings.
def asList(value) {
    if (value == null) {
        return []
    }
    return value.toString().split(',').collect { v -> v.trim() }.findAll { v -> v }
}

// Parse a time parameter into epoch milliseconds (UTC).
// Accepts Unix time in seconds (as used for BigBacter run timestamps) or milliseconds,
// a date (2026-11-09), or a date-time (2026-11-09T14:30:00, optionally with Z or an offset).
def parseTime(value, label) {
    if (value == null || value.toString().trim() == '') {
        return null
    }
    def s = value.toString().trim()
    try {
        if (s.matches('\\d{9,10}')) {
            return s.toLong() * 1000L
        }
        if (s.matches('\\d{13}')) {
            return s.toLong()
        }
        if (s.matches('\\d{4}-\\d{2}-\\d{2}')) {
            return java.time.LocalDate.parse(s).atStartOfDay(java.time.ZoneOffset.UTC).toInstant().toEpochMilli()
        }
        def iso = s.replace(' ', 'T')
        if (iso.matches('.*(Z|[+-]\\d{2}:?\\d{2})')) {
            return java.time.OffsetDateTime.parse(iso).toInstant().toEpochMilli()
        }
        return java.time.LocalDateTime.parse(iso).toInstant(java.time.ZoneOffset.UTC).toEpochMilli()
    }
    catch (Exception _e) {
        error("Could not read --${label} '${s}'. Use Unix time (e.g. 1794182400), a date (2026-11-09) or a date-time (2026-11-09T14:30:00). Times are UTC unless an offset is given.")
    }
}

// Epoch milliseconds -> 2026-11-09T14:30:00Z
def fmtTime(millis) {
    if (millis == null) {
        return ''
    }
    return java.time.Instant.ofEpochMilli(millis.toString().toLong())
        .truncatedTo(java.time.temporal.ChronoUnit.SECONDS)
        .toString()
}

// Bytes -> human readable size
def fmtSize(bytes) {
    if (bytes == null) {
        return ''
    }
    def units = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    def value = bytes.toString().toDouble()
    if (value < 1024) {
        return "${bytes} B"
    }
    def idx = Math.min(units.size() - 1, (Math.log(value) / Math.log(1024)).intValue())
    return String.format('%.1f %s', value / Math.pow(1024, idx), units[idx])
}

// Name, size, last-modified time and type of one path, from a single metadata call.
def fileInfo(path) {
    def attrs = java.nio.file.Files.readAttributes(path, java.nio.file.attribute.BasicFileAttributes)
    def name = path.getFileName().toString()
    if (name.endsWith('/')) {
        name = name.substring(0, name.length() - 1)
    }
    return [
        path : path,
        name : name,
        size : attrs.isDirectory() ? 0L : attrs.size(),
        mtime: attrs.lastModifiedTime()?.toMillis(),
        dir  : attrs.isDirectory(),
    ]
}

// Same as fileInfo, but returns null when the path does not exist.
def fileInfoOrNull(path) {
    try {
        return fileInfo(path)
    }
    catch (java.nio.file.NoSuchFileException _e) {
        return null
    }
    catch (java.io.FileNotFoundException _e) {
        return null
    }
}

// List the direct children of a directory. Returns an empty list if it does not exist.
def listDir(dir) {
    def out = []
    if (dir == null) {
        return out
    }
    def stream = null
    try {
        stream = java.nio.file.Files.newDirectoryStream(dir)
        stream.each { p -> out << fileInfo(p) }
        stream.close()
    }
    catch (java.nio.file.NoSuchFileException _e) {
        return []
    }
    catch (java.nio.file.NotDirectoryException _e) {
        return []
    }
    catch (Exception e) {
        if (stream != null) {
            stream.close()
        }
        throw e
    }
    return out
}

// List every file below a directory (any depth). Returns an empty list if it does not exist.
def listRecursive(dir) {
    def out = []
    listDir(dir).each { info ->
        if (info.dir) {
            out.addAll(listRecursive(info.path))
        }
        else {
            out << info
        }
    }
    return out
}

// Run a closure over a list of items using a small thread pool. Results keep the input order.
def runParallel(items, threads, task) {
    if (!items) {
        return []
    }
    def n = Math.max(1, Math.min(threads.toString().toInteger(), items.size()))
    def pool = java.util.concurrent.Executors.newFixedThreadPool(n)
    try {
        def futures = items.collect { item -> pool.submit({ -> task.call(item) } as java.util.concurrent.Callable) }
        def results = futures.collect { f -> f.get() }
        pool.shutdown()
        return results
    }
    catch (Exception e) {
        pool.shutdownNow()
        throw e
    }
}

// Read a small JSON file (gzipped or not) by streaming it. Used only for config.json and ref.json(.gz).
def readJson(path) {
    def raw = java.nio.file.Files.newInputStream(path)
    try {
        def input = path.getFileName().toString().endsWith('.gz') ? new java.util.zip.GZIPInputStream(raw) : raw
        def json = new groovy.json.JsonSlurper().parse(input, 'UTF-8')
        raw.close()
        return json
    }
    catch (Exception e) {
        raw.close()
        throw e
    }
}

// Compare two numbers given as numbers or strings, so 0.1, 0.10 and 1e-1 are equal.
def sameNumber(a, b) {
    try {
        return new BigDecimal(a.toString()).compareTo(new BigDecimal(b.toString())) == 0
    }
    catch (NumberFormatException _e) {
        return false
    }
}

// Write a tab-separated file from a list of maps.
def writeTsv(path, columns, rows) {
    def lines = [columns.join('\t')]
    rows.each { row ->
        lines << columns.collect { c -> row[c] == null ? '' : row[c].toString().replaceAll('[\t\r\n]', ' ') }.join('\t')
    }
    if (path.toUri().getScheme() == 'file') {
        java.nio.file.Files.createDirectories(path.getParent())
    }
    path.text = lines.join('\n') + '\n'
    return path
}

// Format a list of maps as an aligned text table for the console.
def textTable(columns, rows) {
    def widths = columns.collect { c ->
        ([c.toString().length()] + rows.collect { row -> (row[c] == null ? '' : row[c].toString()).length() }).max()
    }
    def fmt = { values ->
        def cells = []
        values.eachWithIndex { v, i -> cells << v.toString().padRight(widths[i]) }
        return '  ' + cells.join('  ').replaceAll('\\s+$', '')
    }
    def lines = [fmt.call(columns), fmt.call(widths.collect { w -> '-' * w })]
    rows.each { row -> lines << fmt.call(columns.collect { c -> row[c] == null ? '' : row[c] }) }
    return lines.join('\n')
}
