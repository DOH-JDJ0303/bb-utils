#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    bigbacter-db
    Companion tool for the BigBacter pipeline: summarise a BigBacter database, remove
    samples from it, or push the files of a BigBacter run into it.

    There are no processes. Everything runs in the Nextflow head job through the Nextflow
    file API, so no database file is ever staged into a work directory.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { asBool     } from './functions/utils'
include { runSummary } from './functions/summary'
include { runRemove  } from './functions/remove'
include { runPush    } from './functions/push'

def helpMessage() {
    return '''
    bigbacter-db: manage a BigBacter database

    Usage:
      nextflow run main.nf --db <database> --mode <summary|remove|push> [options]

    Required:
      --db <path>              BigBacter database (local path, s3://, gs:// or az://)
      --mode <mode>            summary, remove or push

    General:
      --outdir <path>          Where reports and manifests are written   [bigbacter_db_results]
      --dry_run                remove/push: report what would happen and change nothing
      --threads <n>            Parallel list/copy/delete operations      [16]

    summary:
      --taxa <a,b>             Limit to these taxa (default: all)
      --cluster <a,b>          Limit to these clusters

    remove (always one taxon; choose a selection):
      --taxa <name>            Taxon to remove from (required)
      --cluster <name>         Limit to one cluster. On its own: remove the whole cluster
      --samples <a,b|file>     Sample names, or a file with one per line (needs --cluster)
      --added_after <time>     Samples added at or after this time
      --added_before <time>    Samples added before this time
      --all                    Remove the whole taxon
      --prune_empty <bool>     Delete clusters left with no samples      [true]
                               Times: Unix seconds, 2026-11-09 or 2026-11-09T14:30:00 (UTC)

    push:
      --run <path>             One BigBacter run directory: <outdir>/<timestamp> (required)
      --taxa <a,b>             Limit to these taxa (default: all in the run)
      --overwrite              Replace sample files that are already in the database
      --allow_missing_sig      Push samples that have no signature in the run
    '''.stripIndent()
}

workflow {

    if (asBool(params.help)) {
        log.info(helpMessage())
        exit(0)
    }

    def modes = ['summary', 'remove', 'push']
    if (!params.db) {
        error("--db is required.\n${helpMessage()}")
    }
    if (!(params.mode in modes)) {
        error("--mode must be one of: ${modes.join(', ')}.\n${helpMessage()}")
    }

    def db = file(params.db.toString(), checkIfExists: false)
    def push = params.mode == 'push'
    if (!push && (!db.exists() || !db.isDirectory())) {
        error("--db ${db.toUriString()} does not exist or is not a directory.")
    }

    // Each invocation writes to its own folder, named by Unix time like BigBacter's run folders
    def stamp = java.time.Instant.now().getEpochSecond()
    def outdir = file(params.outdir.toString()).resolve("${stamp}-${params.mode}".toString())

    def opts = [
        db               : db,
        outdir           : outdir,
        threads          : params.threads,
        dry_run          : asBool(params.dry_run),
        taxa             : params.taxa,
        cluster          : params.cluster,
        samples          : params.samples,
        added_after      : params.added_after,
        added_before     : params.added_before,
        all              : asBool(params.all),
        prune_empty      : asBool(params.prune_empty),
        run              : params.run ? file(params.run.toString()) : null,
        overwrite        : asBool(params.overwrite),
        allow_missing_sig: asBool(params.allow_missing_sig),
    ]

    if (params.mode == 'summary') {
        runSummary(opts)
    }
    else if (params.mode == 'remove') {
        runRemove(opts)
    }
    else {
        runPush(opts)
    }
}
