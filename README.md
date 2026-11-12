# bb-utils

A companion to the [BigBacter](https://github.com/nw-page/bigbacter) pipeline for working with a
BigBacter database. It has three modes:

| Mode | What it does |
|---|---|
| `summary` | Describes the database and checks it for problems. Read-only. |
| `remove` | Removes samples from one taxon: by name, by cluster, or by when they were added. |
| `push` | Copies the database files of one BigBacter run into the database. |

It is a Nextflow pipeline with **no processes**. Everything runs in the head job through
Nextflow's file API, so it works on local paths, `s3://`, `gs://` and `az://`, and no database
file is ever staged into a work directory. Listing and deleting only touch file metadata, and
copies within one cloud provider are done by the provider. No container engine is needed.

Requires Nextflow 25.04 or later. Tested on 25.04.8, 25.10.0 and 26.10.5 (strict syntax).

## Quick start

```bash
# What is in the database?
nextflow run DOH-JDJ0303/bb-utils -r main --db s3://bucket/bigbacter_db --mode summary

# Preview a removal, then do it
nextflow run DOH-JDJ0303/bb-utils -r main --db s3://bucket/bigbacter_db --mode remove \
    --taxa Escherichia_coli --cluster 12 --samples SampleA,SampleB --dry_run
nextflow run DOH-JDJ0303/bb-utils -r main --db s3://bucket/bigbacter_db --mode remove \
    --taxa Escherichia_coli --cluster 12 --samples SampleA,SampleB

# Push a finished BigBacter run
nextflow run DOH-JDJ0303/bb-utils -r main --db s3://bucket/bigbacter_db --mode push \
    --run s3://bucket/bigbacter_results/1794182400
```

`remove` and `push` **act immediately**. Add `--dry_run` to see what would happen without
changing anything. Either way, a manifest of every affected file is written before any file is
deleted or copied.

Each invocation writes to `<outdir>/<unix time>-<mode>/` (default `bb_utils_results/`).

## Taxa names

`--taxa` is cleaned the same way BigBacter cleans samplesheet and GAMBIT names: every character
that is not a letter or digit becomes `_`. So `--taxa "Escherichia coli"` and
`--taxa Escherichia_coli` are the same. Case is kept, as in BigBacter, so `escherichia_coli` is a
different taxon.

## summary

```bash
nextflow run DOH-JDJ0303/bb-utils -r main --db <db> --mode summary [--taxa a,b] [--cluster 1,2]
```

| File | Contents |
|---|---|
| `summary_taxa.tsv` | Per taxon: clusters, samples, signatures, size, first and last added, `config.json` values, issue count |
| `summary_clusters.tsv` | Per cluster: samples, reference name, length and contigs, size, first and last added |
| `summary_samples.tsv` | Per sample: taxon, cluster, added time, which files it has, size |
| `summary_issues.tsv` | Problems found (see below) |

Issues reported:

- taxon has no `config.json`, or it has no `clust_dist` (BigBacter v2.2.0+ will not run)
- cluster has no `ref.fa.gz` or no `ref.json.gz`
- cluster has no samples (Floc could reuse its name)
- sample is missing its variant tarball, assembly or signature
- signature with no files in any cluster
- sample present in more than one cluster

A missing signature is expected for samples whose cluster was set in the samplesheet, because
BigBacter does not run Floc for those.

The only file contents read are `config.json` and each cluster's `ref.json.gz`.

## remove

Removals are always limited to one taxon. Choose one selection:

| Selection | Removes |
|---|---|
| `--taxa T --cluster C --samples a,b,c` | Named samples in that cluster. `--samples` can also be a file with one name per line. |
| `--taxa T --cluster C` | The whole cluster. |
| `--taxa T --added_after X --added_before Y` | Samples added in that window. Either bound can be left out. Add `--cluster C` to limit it to one cluster. |
| `--taxa T --all` | The whole taxon, including `config.json`. |

`--samples` and the time bounds can be combined; a sample must then match both.

For each removed sample, its signature, assembly and variant tarball are deleted together. The
signature goes first, so an interrupted run never leaves a sample that Floc can assign to but
that has no cluster files.

**When a sample was added** is the last-modified time of its variant tarball
(`clusters/<cluster>/var/<sample>.tar.gz`), falling back to its assembly. The signature is not
used: Floc rewrites every signature in a taxon on each BigBacter push, so signature times show
the latest push to the taxon, not when the sample arrived.

Times can be Unix seconds (like BigBacter run timestamps), `2026-11-09`, or
`2026-11-09T14:30:00`. They are UTC unless an offset is given. `--added_after` is inclusive and
`--added_before` is exclusive. BigBacter writes to the database at the end of a run, so to remove
what one run added, use that run's timestamp as `--added_after` and the next run's as
`--added_before`.

**Empty clusters are deleted by default.** When a removal leaves a cluster with no samples, its
`ref.fa.gz`, `ref.json.gz` and folder are deleted too. Floc names a new cluster as the highest
existing number plus one, taken from the signatures, and BigBacter treats any existing cluster
folder as an existing cluster. A folder left behind with only a reference could therefore be
picked up by unrelated samples, which would be called against the wrong reference.
`--prune_empty false` keeps the folder and warns.

If a removed sample is the one the cluster's reference was made from, the reference is kept (it
is a separate copy) and a note is printed.

Output: `remove_manifest.tsv`, one row per file with its size, last-modified time and status
(`deleted`, `dry_run` or `failed: ...`). The run exits non-zero if any delete fails.

## push

```bash
nextflow run DOH-JDJ0303/bb-utils -r main --db <db> --mode push --run <outdir>/<timestamp> [--taxa a,b]
```

`--run` is one BigBacter run directory. Files are mapped like this:

| Run directory (`<outdir>/<timestamp>/`) | Database (`<db>/`) |
|---|---|
| `<taxa>/config.json` | `<taxa>/config.json` |
| `<taxa>/sig/<sample>.sig.gz` | `<taxa>/sig/<sample>.sig.gz` |
| `<taxa>/<cluster>/asm/ref.fa.gz` | `<taxa>/clusters/<cluster>/asm/ref.fa.gz` |
| `<taxa>/<cluster>/aux/ref.json.gz` | `<taxa>/clusters/<cluster>/aux/ref.json.gz` |
| `<taxa>/<cluster>/var/<sample>.tar.gz` | `<taxa>/clusters/<cluster>/var/<sample>.tar.gz` |
| `<taxa>/<cluster>/asm/<sample>.fa.gz` | `<taxa>/clusters/<cluster>/asm/<sample>.fa.gz` |

The samples to push are the variant tarballs in the run. BigBacter publishes a tarball only for
samples analysed in that run, so each tarball's path gives the sample, cluster and taxon and no
summary file is needed.

Rules:

- **`config.json`** A taxon is skipped if the run and database `clust_dist` differ, if the
  database has the taxon but no `config.json`, or if the taxon is new and the run has no
  `config.json`. A different `clust_assign_hashfreq` only warns. An existing `config.json` is
  never replaced.
- **References** are never replaced. A cluster is skipped if it is new to the database and the
  run has no `ref.fa.gz`, or if the run and database both have one and they differ in size. The
  second case means the cluster name was reused since the run was made.
- **Samples** are skipped if they are already in the database under a different cluster, or if
  their assembly or signature is missing from the run. `--allow_missing_sig` pushes samples
  without a signature, for runs where clusters were set in the samplesheet.
- **Existing sample files** are left alone unless `--overwrite` is given.
- **Order** `config.json` and references are copied first, then variant tarballs and
  assemblies, then signatures. A sample only becomes visible to Floc once its other files exist.

Output: `push_manifest.tsv`, one row per file with source, destination and status. The run exits
non-zero if any copy fails or if any taxon, cluster or sample was skipped for a problem, so
skipped items are not missed in automated use. Pushing the same run twice copies nothing the
second time.

## BigBacter version needed for push

Push needs a BigBacter version that saves its database files to the run directory. BigBacter
v2.2.0 does not: its run directory only holds the variant tarballs and, for a new taxon,
`config.json`. Push also needs:

| File | Publish to | Note |
|---|---|---|
| Signatures (`FLOC` `sigs/*.sig.gz`) | `<timestamp>/<taxa>/sig/` | Floc outputs every signature in the taxon; this tool only pushes those with a variant tarball in the run. |
| Reference (`ref.fa.gz`) | `<timestamp>/<taxa>/<cluster>/asm/` | Best published for every cluster in the run, not only new ones, so push can confirm the run and database used the same reference. |
| Reference metadata (`ref.json.gz` or `ref.json`) | `<timestamp>/<taxa>/<cluster>/aux/` | Either name is accepted and copied as is. |
| Sample assembly (`STD_ASSEMBLY` output) | `<timestamp>/<taxa>/<cluster>/asm/<sample>.fa.gz` | The gzipped copy that goes into the database. It replaces the uncompressed `<sample>/asm/<sample>.fa` that `SHOVILL` and `NCBI_DATASETS` publish today. |

Optional: publishing `config.json` to `<timestamp>/<taxa>/` on every run, not only the run that
creates the taxon, would let push check `clust_dist` for existing taxa too.

## Things to know

- **Deletes are permanent** unless the bucket has versioning. The manifest lists what was
  removed.
- **Last-modified times are reset** when files are rewritten: copying or syncing a database to a
  new location, or pushing with `--overwrite`. After a migration every sample looks like it was
  added on the migration date.
- **No locking.** Do not run `remove` or `push` while a BigBacter `--push true` run is writing to
  the same database.
- **Permissions.** The head job needs list, read, write and delete on the database.
- **Stale runs.** If a run has no `ref.fa.gz` for a cluster, push cannot tell whether the
  database cluster of that name is still the one the run used. This matters only if clusters were
  removed and their names reused between the run and the push.
- **Scale.** Listing is done by the head job with `--threads` (default 16) parallel calls: one
  per taxon plus three per cluster. Very large databases take longer and need more head-job
  memory, but no data is transferred.

## Parameters

| Parameter | Default | Modes | Description |
|---|---|---|---|
| `--db` | | all | BigBacter database |
| `--mode` | | all | `summary`, `remove` or `push` |
| `--outdir` | `bb_utils_results` | all | Reports and manifests |
| `--threads` | `16` | all | Parallel list, copy and delete operations |
| `--dry_run` | `false` | remove, push | Report only, change nothing |
| `--taxa` | | all | remove: one taxon (required). summary, push: optional, comma-separated |
| `--cluster` | | summary, remove | remove: one cluster. summary: optional, comma-separated |
| `--samples` | | remove | Comma-separated names or a file with one per line |
| `--added_after` | | remove | Samples added at or after this time |
| `--added_before` | | remove | Samples added before this time |
| `--all` | `false` | remove | Remove the whole taxon |
| `--prune_empty` | `true` | remove | Delete clusters left with no samples |
| `--run` | | push | BigBacter run directory, `<outdir>/<timestamp>` |
| `--overwrite` | `false` | push | Replace sample files already in the database |
| `--allow_missing_sig` | `false` | push | Push samples with no signature in the run |

## Tests

```bash
tests/run_tests.sh            # uses `nextflow` on the PATH, or pass a binary as the first argument
```

Builds a synthetic local database and run directory (`tests/make_test_data.sh`) and checks all
three modes, including dry runs, refusals and skip rules. Cloud paths use the same code but are
not covered by these tests.
