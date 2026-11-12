#!/usr/bin/env bash
# Builds a small synthetic BigBacter database and run directory for testing bb-utils.
# Usage: tests/make_test_data.sh <target_dir>
set -euo pipefail

out="${1:?usage: make_test_data.sh <target_dir>}"
rm -rf "$out"
db="$out/db"
run="$out/outdir/1762732800"

# mk <path> <mtime> [content]  : write a gzipped file with a fixed last-modified time
mk() {
    mkdir -p "$(dirname "$1")"
    printf '%s\n' "${3:-$(basename "$1")}" | gzip -c > "$1"
    touch -d "$2" "$1"
}

# sample <taxa> <cluster> <sample> <added> <sig_time>
sample() {
    mk "$db/$1/clusters/$2/asm/$3.fa.gz"  "$4" ">$3 ACGTACGTACGT"
    mk "$db/$1/clusters/$2/var/$3.tar.gz" "$4" "snippy $3"
    mk "$db/$1/sig/$3.sig.gz"             "$5" "[{\"name\":\"$3\",\"filename\":\"$2\"}]"
}

# ref <taxa> <cluster> <ref_sample> <time>
ref() {
    mk "$db/$1/clusters/$2/asm/ref.fa.gz"   "$4" ">ref ACGTACGTACGT"
    mk "$db/$1/clusters/$2/aux/ref.json.gz" "$4" "{\"length\": 4840898, \"n_contigs\": 3, \"name\": \"$3\"}"
}

# ---------------------------------------------------------------------------------------
#  Database
# ---------------------------------------------------------------------------------------
# Signatures all carry the time of the latest push (Floc rewrites every signature each run)
SIG="2026-06-01T12:00:00Z"

t=Escherichia_coli
mkdir -p "$db/$t"
cat > "$db/$t/config.json" <<EOF
{
  "clust_dist": 0.1,
  "clust_assign_hashfreq": 0.6,
  "floc_version": "1.2.4",
  "bigbacter_version": "2.2.0"
}
EOF
ref    $t 1 S1 "2026-01-10T12:00:00Z"
sample $t 1 S1 "2026-01-10T12:00:00Z" "$SIG"
sample $t 1 S2 "2026-01-10T12:00:00Z" "$SIG"
sample $t 1 S3 "2026-03-05T12:00:00Z" "$SIG"
ref    $t 2 S4 "2026-01-10T12:00:00Z"
sample $t 2 S4 "2026-01-10T12:00:00Z" "$SIG"
sample $t 2 S5 "2026-03-05T12:00:00Z" "$SIG"
ref    $t 3 S6 "2026-06-01T12:00:00Z"
sample $t 3 S6 "2026-06-01T12:00:00Z" "$SIG"

# A pre-v2.2.0 taxon with problems: no config.json, a sample missing its var, an orphan
# signature, and a cluster with a reference but no samples
t=Klebsiella_pneumoniae
ref    $t 1 K1 "2025-11-01T12:00:00Z"
sample $t 1 K1 "2025-11-01T12:00:00Z" "2025-11-02T12:00:00Z"
sample $t 1 K2 "2025-11-02T12:00:00Z" "2025-11-02T12:00:00Z"
rm "$db/$t/clusters/1/var/K2.tar.gz"
mk "$db/$t/sig/K9.sig.gz" "2025-11-02T12:00:00Z"
ref    $t 2 K7 "2025-11-01T12:00:00Z"

# ---------------------------------------------------------------------------------------
#  BigBacter run directory (<outdir>/<timestamp>), using the agreed paths
# ---------------------------------------------------------------------------------------
NOW="2026-11-09T12:00:00Z"

# rsample <taxa> <cluster> <sample> [nosig] : files BigBacter publishes for a sample it analysed
rsample() {
    mk "$run/$1/$2/asm/$3.fa.gz"     "$NOW" ">$3 ACGTACGTACGT"
    mk "$run/$1/$2/var/$3.tar.gz"    "$NOW" "snippy $3"
    mkdir -p "$run/$3/taxa" && echo "predicted.name" > "$run/$3/taxa/$3_gambit.csv"
    if [ "${4:-}" != "nosig" ]; then
        mk "$run/$1/sig/$3.sig.gz"   "$NOW" "[{\"name\":\"$3\",\"filename\":\"$2\"}]"
    fi
}
# rref <taxa> <cluster> <ref_sample> [content]
rref() {
    mk "$run/$1/$2/asm/ref.fa.gz"   "$NOW" "${4:->ref ACGTACGTACGT}"
    mk "$run/$1/$2/aux/ref.json.gz" "$NOW" "{\"length\": 5000000, \"n_contigs\": 2, \"name\": \"$3\"}"
}

t=Escherichia_coli
mkdir -p "$run/$t/cluster" && echo "sid,cluster" > "$run/$t/cluster/clusters.csv"
rsample $t 1 N1            # new sample, existing cluster
rsample $t 1 S2            # already in the database
rsample $t 4 N2            # new cluster with its reference
rsample $t 4 N3
rref    $t 4 N2
rsample $t 2 N4 nosig      # no signature in the run
rsample $t 3 N7            # run's reference differs from the database's
rref    $t 3 N7 ">a different reference ACGTACGTACGTACGTACGT"
rsample $t 2 S1            # already in the database under cluster 1
mk "$run/$t/1/summary/1762732800-$t-1_summary.csv" "$NOW"

t=Salmonella_enterica      # taxon new to the database
mkdir -p "$run/$t"
cat > "$run/$t/config.json" <<EOF
{"clust_dist": 0.1, "clust_assign_hashfreq": 0.6, "floc_version": "1.2.4", "bigbacter_version": "2.2.0"}
EOF
rsample $t 1 N5
rref    $t 1 N5

t=Staphylococcus_aureus    # taxon new to the database, but the run has no config.json
rsample $t 1 N6
rref    $t 1 N6

mkdir -p "$run/multiqc" && echo "<html>" > "$run/multiqc/multiqc_report.html"
echo "id,run,status" > "$run/1762732800-summary.csv"

echo "Test data written to $out"
