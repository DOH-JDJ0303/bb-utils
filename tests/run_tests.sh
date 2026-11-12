#!/usr/bin/env bash
# End-to-end tests for bigbacter-db against a synthetic local database.
# Usage: tests/run_tests.sh [nextflow_binary]
set -uo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
nxf="${1:-nextflow}"
work="$(mktemp -d)"
cd "$work"
export NXF_ANSI_LOG=false

pass=0
fail=0
ok()   { pass=$((pass + 1)); echo "  ok    $1"; }
bad()  { fail=$((fail + 1)); echo "  FAIL  $1"; }
check() { if eval "$2"; then ok "$1"; else bad "$1"; fi; }

reset() { "$here/make_test_data.sh" "$work/data" > /dev/null; }
run()   { "$nxf" run "$here/../main.nf" --db data/db "$@" > last.log 2>&1; echo $?; }
latest() { ls -d bigbacter_db_results/*-"$1" | sort | tail -1; }
has()   { [ -e "data/db/$1" ]; }
hasnt() { [ ! -e "data/db/$1" ]; }
state() { find data/db -type f | sort | xargs md5sum | md5sum; }

E=Escherichia_coli

echo "summary"
reset
rc=$(run --mode summary)
out=$(latest summary)
check "exits 0"                          '[ "$rc" = 0 ]'
check "2 taxa"                           '[ "$(tail -n +2 $out/summary_taxa.tsv | wc -l)" = 2 ]'
check "5 clusters"                       '[ "$(tail -n +2 $out/summary_clusters.tsv | wc -l)" = 5 ]'
check "8 samples"                        '[ "$(tail -n +2 $out/summary_samples.tsv | wc -l)" = 8 ]'
check "4 issues"                         '[ "$(tail -n +2 $out/summary_issues.tsv | wc -l)" = 4 ]'
check "added time comes from var, not sig" 'grep -q "S1	2026-01-10T12:00:00Z" $out/summary_samples.tsv'
check "reference name is read"           'grep -q "^$E	1	3	S1	4840898	3" $out/summary_clusters.tsv'
rc=$(run --mode summary --taxa "Escherichia coli")
check "taxa name is cleaned like BigBacter" '[ "$rc" = 0 ] && [ "$(tail -n +2 $(latest summary)/summary_taxa.tsv | wc -l)" = 1 ]'

echo "remove: dry run changes nothing"
reset; before=$(state)
rc=$(run --mode remove --taxa $E --cluster 1 --samples S1,S3 --dry_run)
check "exits 0"                          '[ "$rc" = 0 ]'
check "database unchanged"               '[ "$before" = "$(state)" ]'
check "manifest lists 6 files"           '[ "$(grep -c dry_run $(latest remove)/remove_manifest.tsv)" = 6 ]'

echo "remove: by name"
reset
rc=$(run --mode remove --taxa $E --cluster 1 --samples S1,S3)
check "exits 0"                          '[ "$rc" = 0 ]'
check "S1 and S3 files are gone"         'hasnt $E/sig/S1.sig.gz && hasnt $E/clusters/1/var/S1.tar.gz && hasnt $E/clusters/1/asm/S3.fa.gz && hasnt $E/sig/S3.sig.gz'
check "S2 is kept"                       'has $E/sig/S2.sig.gz && has $E/clusters/1/var/S2.tar.gz && has $E/clusters/1/asm/S2.fa.gz'
check "reference is kept"                'has $E/clusters/1/asm/ref.fa.gz && has $E/clusters/1/aux/ref.json.gz'
check "other clusters untouched"         'has $E/clusters/2/var/S4.tar.gz && has $E/clusters/3/var/S6.tar.gz'
check "warns that S1 was the reference"  'grep -q "S1 is the reference" last.log'

echo "remove: sample names from a file"
reset; printf 'S4\nS5\n' > names.txt
rc=$(run --mode remove --taxa $E --cluster 2 --samples names.txt)
check "emptied cluster is deleted"       '[ "$rc" = 0 ] && hasnt $E/clusters/2 && hasnt $E/sig/S4.sig.gz && hasnt $E/sig/S5.sig.gz'

echo "remove: by time"
reset
rc=$(run --mode remove --taxa $E --added_after 2026-03-01 --added_before 2026-07-01)
check "exits 0"                          '[ "$rc" = 0 ]'
check "S3, S5, S6 removed"               'hasnt $E/sig/S3.sig.gz && hasnt $E/clusters/2/var/S5.tar.gz && hasnt $E/sig/S6.sig.gz'
check "earlier samples kept"             'has $E/clusters/1/var/S1.tar.gz && has $E/clusters/2/var/S4.tar.gz && has $E/sig/S2.sig.gz'
check "emptied cluster 3 is deleted"     'hasnt $E/clusters/3'
check "other taxon untouched"            'has Klebsiella_pneumoniae/clusters/1/var/K1.tar.gz'
reset
rc=$(run --mode remove --taxa $E --added_after 1767225600 --cluster 1)
check "Unix time and --cluster limit"    '[ "$rc" = 0 ] && hasnt $E/clusters/1/var/S1.tar.gz && has $E/clusters/2/var/S4.tar.gz'

echo "remove: --prune_empty false keeps the reference"
reset
rc=$(run --mode remove --taxa $E --cluster 3 --prune_empty false)
check "samples gone, reference kept"     '[ "$rc" = 0 ] && hasnt $E/clusters/3/var/S6.tar.gz && has $E/clusters/3/asm/ref.fa.gz'

echo "remove: whole cluster and whole taxon"
reset
rc=$(run --mode remove --taxa $E --cluster 1)
check "cluster 1 deleted"                '[ "$rc" = 0 ] && hasnt $E/clusters/1 && hasnt $E/sig/S1.sig.gz && has $E/sig/S4.sig.gz && has $E/config.json'
rc=$(run --mode remove --taxa $E --all)
check "taxon deleted"                    '[ "$rc" = 0 ] && hasnt $E && has Klebsiella_pneumoniae'

echo "remove: refuses unsafe or unclear requests"
reset; before=$(state)
rc=$(run --mode remove);                                   check "no --taxa"            '[ "$rc" != 0 ]'
rc=$(run --mode remove --taxa $E);                         check "no selection"         '[ "$rc" != 0 ]'
rc=$(run --mode remove --taxa $E --samples S1);            check "--samples without --cluster" '[ "$rc" != 0 ]'
rc=$(run --mode remove --taxa escherichia_coli --all);     check "wrong-case taxon"     '[ "$rc" != 0 ]'
rc=$(run --mode remove --taxa $E --cluster 99);            check "unknown cluster"      '[ "$rc" != 0 ]'
rc=$(run --mode remove --taxa $E --added_after yesterday); check "unreadable time"      '[ "$rc" != 0 ]'
rc=$(run --mode remove --taxa $E --added_after 2026-07-01 --added_before 2026-01-01); check "reversed window" '[ "$rc" != 0 ]'
check "database unchanged"               '[ "$before" = "$(state)" ]'

echo "push: dry run changes nothing"
reset; before=$(state)
rc=$(run --mode push --run data/outdir/1762732800 --dry_run)
check "exits 0"                          '[ "$rc" = 0 ]'
check "database unchanged"               '[ "$before" = "$(state)" ]'

echo "push"
reset; refsum=$(md5sum < data/db/$E/clusters/3/asm/ref.fa.gz); s2=$(stat -c %Y data/db/$E/clusters/1/var/S2.tar.gz)
rc=$(run --mode push --run data/outdir/1762732800)
check "exits non-zero because items were skipped" '[ "$rc" != 0 ]'
check "new sample in existing cluster"   'has $E/clusters/1/var/N1.tar.gz && has $E/clusters/1/asm/N1.fa.gz && has $E/sig/N1.sig.gz'
check "new cluster with its reference"   'has $E/clusters/4/asm/ref.fa.gz && has $E/clusters/4/aux/ref.json.gz && has $E/clusters/4/var/N3.tar.gz && has $E/sig/N2.sig.gz'
check "new taxon with config.json"       'has Salmonella_enterica/config.json && has Salmonella_enterica/clusters/1/asm/ref.fa.gz && has Salmonella_enterica/sig/N5.sig.gz'
check "existing sample not rewritten"    '[ "$s2" = "$(stat -c %Y data/db/$E/clusters/1/var/S2.tar.gz)" ]'
check "sample without signature skipped" 'hasnt $E/clusters/2/var/N4.tar.gz && hasnt $E/clusters/2/asm/N4.fa.gz'
check "sample in another cluster skipped" 'hasnt $E/clusters/2/var/S1.tar.gz'
check "mismatched reference: cluster skipped, reference untouched" 'hasnt $E/clusters/3/var/N7.tar.gz && hasnt $E/sig/N7.sig.gz && [ "$refsum" = "$(md5sum < data/db/$E/clusters/3/asm/ref.fa.gz)" ]'
check "taxon without config.json skipped" 'hasnt Staphylococcus_aureus'
check "non-database files not copied"    'hasnt $E/cluster && hasnt $E/clusters/1/summary && hasnt multiqc'
before=$(state)
rc=$(run --mode push --run data/outdir/1762732800)
check "second push copies nothing"       '[ "$before" = "$(state)" ] && grep -q "^0 file(s) to copy" last.log'
rc=$(run --mode push --run data/outdir/1762732800 --taxa "Salmonella enterica" --overwrite)
check "--taxa and --overwrite"           '[ "$rc" = 0 ] && grep -q "^3 file(s) to copy" last.log'
rc=$(run --mode push --run data/outdir/1762732800 --allow_missing_sig)
check "--allow_missing_sig"              'has $E/clusters/2/var/N4.tar.gz && has $E/clusters/2/asm/N4.fa.gz && hasnt $E/sig/N4.sig.gz'

echo "push: config.json must match"
reset; sed -i 's/"clust_dist": 0.1/"clust_dist": 0.03/' data/db/$E/config.json; before=$(state)
rc=$(run --mode push --run data/outdir/1762732800 --taxa $E)
check "existing taxon, run has no config: allowed" 'has $E/sig/N1.sig.gz'
reset; printf '{"clust_dist": 0.03}\n' > data/outdir/1762732800/$E/config.json; before=$(state)
rc=$(run --mode push --run data/outdir/1762732800 --taxa $E)
check "different clust_dist: taxon skipped" '[ "$rc" != 0 ] && [ "$before" = "$(state)" ]'
reset
rc=$(run --mode push --run data/outdir);  check "wrong --run level is refused" '[ "$rc" != 0 ]'

echo
echo "$pass passed, $fail failed"
cd / && rm -rf "$work"
[ "$fail" = 0 ]
