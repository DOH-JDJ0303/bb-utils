# bb-utils

A companion workflow to [BigBacter](https://github.com/nw-page/bigbacter) for managing a
BigBacter database. It has three modes:

| Mode | What it does |
|---|---|
| `summary` | Creates a summary of the database (read-only). |
| `remove` | Removes samples by name, cluster, or date added. |
| `push` | Pushes files from a BigBacter run to the database. |

bb-utils requires Nextflow 25.04 or later and BigBacter v2.2.1 or later.

## Usage

```bash
nextflow run nw-page/bb-utils \
    -r main \
    --db s3://bucket/bigbacter_db \
    --mode summary
```

Add `--dry_run` to `remove` and `push` to test a command before executing it.

## Documentation

Full usage, parameters and outputs are described on the
[bb-utils page](https://nw-page.github.io/bigbacter-docs/docs/auxillary/pages/bb-utils/) of the
BigBacter docs.