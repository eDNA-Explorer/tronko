#!/usr/bin/env bash
# --hits tied on the multi-tree fixture, whose three tree copies hold the same sequences (copy 0
# unsuffixed, copy 1 every second leaf with _1, copy 2 with _2), so a read ties leaves in two or
# three trees:
#   - each case reproduces its golden (goldens-hits-tied, recorded at one thread) at every thread
#     count, and --hits primary reproduces the default golden;
#   - BWA-MEM3 gets exactly -z 1 -h 2147483647 -u after the pinned options, and they change no
#     SAM field Tronko reads (QNAME FLAG RNAME POS CIGAR RNEXT) but the XA text;
#   - every read's candidate trees (from --trace-read) are, in ascending order, its --hits
#     primary trees and the trees of its ties, recomputed here from the aligner's own SAM of the
#     batch: per mate the primary record's leaf and each XA entry with the record's AS on the
#     strand Tronko aligns (rule H8); for a pair the trees both mates tie on, if any;
#   - it is refused off the node-store path and for an unknown mode.
# Environment: golden_lib.sh (TRONKO_TESTS_RECORD=1 records the goldens).
set -uo pipefail
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/golden_lib.sh"

MT=${TRONKO_MULTITREE_DIR:-"$TA/tests/data/multitree"}
DATA="$REPO_ROOT/tests/data/assignment"
for f in reference_tree.trkb multitree.fasta single_F.fasta single_R.fasta; do
	[[ -f "$MT/$f" ]] || { echo "SKIP: missing fixture file $MT/$f"; exit 77; }
done
for f in "$DATA/paired_2000_1.fasta" "$DATA/paired_2000_2.fasta"; do
	[[ -f $f ]] || { echo "SKIP: missing fixture file $f"; exit 77; }
done
# Recorded from a cap-10 build only.
[[ $MATCH_CAP == 10 ]] || { echo "SKIP: --hits tied goldens exist for cap 10 only"; exit 77; }
mem3_index "$MT/multitree.fasta" || exit 2

BASE=(-r -f "$MT/reference_tree.trkb" -a "$MT/multitree.fasta" -w -6 --Cinterval 10)
GOLDENS="$MT/goldens-hits-tied"
DEFAULTS=$(golden_dir "$MT/goldens")
MEM=(--compat=bwa-mem -K 100000000 --dedup-reads off --dedup off --ks-dedup off --cohort-slices 0 -v 1)
TIED=(-z 1 -h 2147483647 -u)
mkdir -p "$GOLDENS"
echo "hits tied: multi-tree fixture, threads $CORES$([[ $RECORD == 1 ]] && echo ', recording goldens')"

# A stand-in bwa-mem3 that logs each command line, then runs the real one.
REAL=$TRONKO_BWA_MEM3
mkdir -p "$TMP_DIR/bin"
cat >"$TMP_DIR/bin/bwa-mem3" <<EOF
#!/usr/bin/env bash
printf '%s\n' "\$*" >> "$TMP_DIR/argv.log"
exec "$REAL" "\$@"
EOF
chmod +x "$TMP_DIR/bin/bwa-mem3"

# batch_fasta <out> <fasta> <reverse 0|1>: the reads as Tronko hands them to the aligner (named
# by their index, upper case, reverse-complemented for -v and the second mate of -z), and
# <out>.names: index, read name
batch_fasta() {
	awk -v rc="$3" -v names="$1.names" '
		function flush() { if (seq == "") return; s = toupper(seq)
			if (rc) { r = ""; for (i = length(s); i > 0; i--) { c = substr(s, i, 1)
				r = r (c == "A" ? "T" : c == "C" ? "G" : c == "G" ? "C" : c == "T" ? "A" : c) }; s = r }
			print n "\t" name > names; print ">" n++; print s; seq = "" }
		BEGIN { n = 0 }
		/^>/ { flush(); name = substr($1, 2); next } { seq = seq $0 } END { flush() }' "$2" > "$1"
}

# expected_trees <sam> <names> <primary trace> <paired 0|1>: "read<TAB>trees" for every read
expected_trees() {
	awk -F'\t' -v paired="$4" '
		function bit(x, b) { return int(x / b) % 2 }
		function tree(name) { return name ~ /_1$/ ? 1 : name ~ /_2$/ ? 2 : 0 }
		function mark(t, m) { if (!bit(from[t], 2 ^ m)) from[t] += 2 ^ m }
		FILENAME == ARGV[1] {
			q = $1; f = $2; s = paired ? 2 * q + bit(f, 128) : q
			p10[s] += bit(f, 16)
			p20[s] += bit(f, 32)
			if (!bit(f, 256) && !bit(f, 2048) && !bit(f, 4)) {
				rname[s] = $3; as[s] = ""; xa[s] = ""
				for (i = 12; i <= NF; i++) {
					if ($i ~ /^AS:i:/) as[s] = substr($i, 6) + 0
					if ($i ~ /^XA:Z:/) xa[s] = substr($i, 6)
				}
			}
			next
		}
		FILENAME == ARGV[2] { name[$1] = $2; if ($1 + 1 > nreads) nreads = $1 + 1; next }
		$1 == "trace" { prim[$2] = $4 }
		END {
			for (q = 0; q < nreads; q++) {
				delete from
				for (m = 0; m < (paired ? 2 : 1); m++) {
					s = paired ? 2 * q + m : q
					rc = (paired ? (m == 0 ? p10[2 * q] : p20[2 * q]) : p10[s]) % 2
					if (!(s in rname)) continue
					mark(tree(rname[s]), m)
					n = split(xa[s], e, ";")
					for (k = 1; k <= n; k++) {
						if (e[k] == "") continue
						split(e[k], g, ",")
						if (g[5] + 0 != as[s]) continue
						if ((substr(g[2], 1, 1) == "-") != rc) continue
						mark(tree(g[1]), m)
					}
				}
				both = 0
				for (t in from) if (from[t] == 3) both = 1
				delete want
				for (t in from) if (!both || from[t] == 3) want[t] = 1
				p = prim[name[q]]
				if (p != "-" && p != "") { np = split(p, pt, ","); for (k = 1; k <= np; k++) want[pt[k]] = 1 }
				out = ""
				for (t = 0; t <= 2; t++) if (t in want) out = out (out == "" ? "" : ",") t
				print name[q] "\t" (out == "" ? "-" : out)
			}
		}' "$1" "$2" "$3"
}

run_case() {
	local c=$1 paired=$2 n out trace bad
	shift 2
	# --hits primary, spelled out, is the default
	if ! TRONKO_BWA_MEM3=$TMP_DIR/bin/bwa-mem3 run_assign "$TMP_DIR/$c-primary.log" "${BASE[@]}" --hits primary \
		--number-of-cores 4 "$@" --trace-read '*' --trace-file "$TMP_DIR/$c-primary.trace" -o "$TMP_DIR/$c-primary.tsv"; then
		FAILS=$((FAILS + 1))
		return
	fi
	[[ $RECORD == 1 ]] || check_output "$c" "$DEFAULTS/expected_$c.tsv" "$TMP_DIR/$c-primary.tsv" "$c, --hits primary"
	for n in $CORES; do
		out=$TMP_DIR/$c-c$n.tsv
		trace=$TMP_DIR/$c-c$n.trace
		if ! TRONKO_BWA_MEM3=$TMP_DIR/bin/bwa-mem3 run_assign "$TMP_DIR/$c-c$n.log" "${BASE[@]}" --hits tied \
			--number-of-cores "$n" "$@" --trace-read '*' --trace-file "$trace" -o "$out"; then
			FAILS=$((FAILS + 1))
			continue
		fi
		check_output "$c" "$GOLDENS/expected_$c.tsv" "$out" "$c, --hits tied, --number-of-cores $n"
	done
	# the aligner's own SAM of the batch, with the tied options and without
	local fa=("$TMP_DIR/$c.1.fa")
	if [[ $paired == 1 ]]; then
		batch_fasta "$TMP_DIR/$c.1.fa" "$DATA/paired_2000_1.fasta" 0
		batch_fasta "$TMP_DIR/$c.2.fa" "$DATA/paired_2000_2.fasta" 1
		fa+=("$TMP_DIR/$c.2.fa")
	else
		local src=${*: -1} rc=0
		[[ " $* " == *" -v "* ]] && rc=1
		batch_fasta "$TMP_DIR/$c.1.fa" "$src" "$rc"
	fi
	"$REAL" mem "${MEM[@]}" "${TIED[@]}" -t 4 "$MT/multitree.fasta" "${fa[@]}" 2>/dev/null | grep -v '^@' > "$TMP_DIR/$c.tied.sam"
	"$REAL" mem "${MEM[@]}" -t 4 "$MT/multitree.fasta" "${fa[@]}" 2>/dev/null | grep -v '^@' > "$TMP_DIR/$c.sam"
	if ! cmp -s <(cut -f1-4,6,7 "$TMP_DIR/$c.sam") <(cut -f1-4,6,7 "$TMP_DIR/$c.tied.sam"); then
		echo "FAIL: $c: the tied options changed SAM fields Tronko reads" >&2
		FAILS=$((FAILS + 1))
	fi
	expected_trees "$TMP_DIR/$c.tied.sam" "$TMP_DIR/$c.1.fa.names" "$TMP_DIR/$c-primary.trace" "$paired" | sort > "$TMP_DIR/$c.want"
	awk -F'\t' '$1 == "trace" { print $2 "\t" $4 }' "$trace" | sort > "$TMP_DIR/$c.got"
	bad=$(comm -3 "$TMP_DIR/$c.want" "$TMP_DIR/$c.got" | wc -l)
	local gained
	gained=$(paste <(sort "$TMP_DIR/$c.want") <(awk -F'\t' '$1 == "trace" { print $2 "\t" $4 }' "$TMP_DIR/$c-primary.trace" | sort) \
		| awk -F'\t' '{ n2 = $2 == "-" ? 0 : split($2, x, ","); n4 = $4 == "-" ? 0 : split($4, y, ",") } n2 > n4' | wc -l)
	echo "  $c: $gained reads gain trees by their ties"
	if [[ $bad -ne 0 || $gained -eq 0 ]]; then
		echo "FAIL: $c: $bad lines differ between the recomputed and the traced candidate trees" >&2
		comm -3 "$TMP_DIR/$c.want" "$TMP_DIR/$c.got" | head -5 >&2
		FAILS=$((FAILS + 1))
	fi
}

run_case mt_single_F 0 -s -g "$MT/single_F.fasta"
run_case mt_unpaired_R 0 -s -v -g "$MT/single_R.fasta"
run_case mt_paired 1 -p -z -1 "$DATA/paired_2000_1.fasta" -2 "$DATA/paired_2000_2.fasta"

# Batch size: BWA-MEM3 picks the primary among equal hits by the read's place in the batch, so
# --hits primary changes with -L; the set of tied hits does not, and a plain (not chimeric)
# single-end read gets the same row at -L 400 as in one batch. (A chimera's tied set is that of
# whichever segment is primary, so it still moves.)
if [[ $RECORD != 1 ]]; then
	for h in primary tied; do
		run_assign "$TMP_DIR/L400-$h.log" "${BASE[@]}" --hits "$h" -L 400 --number-of-cores 4 -s -g "$MT/single_F.fasta" \
			-o "$TMP_DIR/L400-$h.tsv" || FAILS=$((FAILS + 1))
	done
	plain_moved() {   # plain_moved <one batch> <-L 400>: plain reads whose row differs
		paste "$1" "$2" | awk -F'\t' 'NR > 1 && $1 ~ /_plain_/ && $0 != "" { half = NF / 2; for (i = 1; i <= half; i++) if ($i != $(i + half)) { n++; break } } END { print n + 0 }'
	}
	moved_primary=$(plain_moved "$DEFAULTS/expected_mt_single_F.tsv" "$TMP_DIR/L400-primary.tsv")
	moved_tied=$(plain_moved "$GOLDENS/expected_mt_single_F.tsv" "$TMP_DIR/L400-tied.tsv")
	echo "  -L 400: plain reads whose row moves: $moved_primary with --hits primary, $moved_tied with --hits tied"
	if [[ $moved_tied -ne 0 || $moved_primary -eq 0 ]]; then
		echo "FAIL: -L 400 moved $moved_tied plain reads under --hits tied ($moved_primary under primary)" >&2
		FAILS=$((FAILS + 1))
	fi
fi

# The command lines: the pinned options, then the tied ones exactly when --hits tied
pinned="mem ${MEM[*]}"
plain=$(grep -c "^$pinned -t " "$TMP_DIR/argv.log")
tied=$(grep -c "^$pinned ${TIED[*]} -t " "$TMP_DIR/argv.log")
total=$(grep -c '^mem ' "$TMP_DIR/argv.log")
if [[ $plain -eq 0 || $tied -eq 0 || $((plain + tied)) -ne $total ]]; then
	echo "FAIL: aligner command lines: $plain plain, $tied with the tied options, $total in all" >&2
	FAILS=$((FAILS + 1))
fi

for env in "TRONKO_NODESTORE=0" "TRONKO_NODESTORE=1"; do
	extra=()
	[[ $env == TRONKO_NODESTORE=1 ]] && extra=(--enable-pruning)
	if env "$env" "$ASSIGN" "${BASE[@]}" --hits tied "${extra[@]}" -s -g "$MT/single_F.fasta" -o "$TMP_DIR/refused.tsv" \
		>/dev/null 2>&1; then
		echo "FAIL: --hits tied off the node-store path ($env ${extra[*]}) was not refused" >&2
		FAILS=$((FAILS + 1))
	fi
done
if "$ASSIGN" "${BASE[@]}" --hits all -s -g "$MT/single_F.fasta" -o "$TMP_DIR/refused.tsv" >/dev/null 2>&1; then
	echo "FAIL: --hits all was not refused" >&2
	FAILS=$((FAILS + 1))
fi

((FAILS == 0)) || { echo "hits tied: $FAILS failure(s)" >&2; exit 1; }
echo "hits tied: PASS"
