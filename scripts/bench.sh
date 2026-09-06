#!/usr/bin/env bash
# Reproduces the tables in docs/REPORT.md.  Usage: scripts/bench.sh [nefvol-binary]
set -euo pipefail
cd "$(dirname "$0")/.."
BIN=${1:-./target/release/nefvol}
stat() { "$BIN" volume --input "$1" --threads "$2" --stats-json 2>&1 >/dev/null | python3 -c '
import json,sys; j=json.load(sys.stdin)
print("%-16s H=%4d facets=%5d visited=%8d emitted=%7d bases=%8d terms=%6d depth=%4d stackdirs=%5d knots=%6d prep=%5dms search=%7dms rss=%6dkB" % (
 sys.argv[1], j["hyperplanes"], j["facet_tasks"], j["vertices_visited"], j["vertices_emitted"], j["bases"], j["terms"], j["max_depth"], j["max_stack_dirs"], j["knots"], j["prep_ms"], j["search_ms"], j["peak_rss_kb"]))' "$(basename "$1" .ine)"; }
echo "# single thread, growing instances (memory flatness)"
for f in boxes3d_20 boxes3d_40 boxes3d_80 boxes3d_160 boxes3d_320; do stat inputs/$f.ine 1; done
echo "# other families, single thread"
for f in rand3d_20 rand3d_60 rand4d_8 rand4d_24 boxes4d_12 boxes5d_10 boxes6d_6 cones3d_30x3 grid3d_4; do stat inputs/$f.ine 1; done
echo "# thread scaling on boxes3d_160"
for t in 1 2 4 8 16; do echo -n "threads=$t  "; stat inputs/boxes3d_160.ine $t; done
echo "# subtree histogram (complete subtrees, budget 0) on rand3d_20"
"$BIN" volume --input inputs/rand3d_20.ine --budget 0 --threads 1 --stats-json 2>&1 >/dev/null | python3 -c 'import json,sys; j=json.load(sys.stdin); print(j["subtree_histogram"])'
