# Real-reference tests of tronko-assign (../README.md, "Real-reference goldens"). The read pairs and
# goldens are committed under ../data/real/; the reference files are fetched by URL and checked
# against ../manifests/<marker>.sha256.
#   make test-real         this checkout on production's real references, against the committed
#                          one-thread goldens (test-real.sh)
#   make verify-branches   test-real on each of BRANCHES, one summary table (verify-branches.sh)
#   make goldens           make the goldens again from GOLDEN_COMMIT at one thread and
#                          compare them with the committed ones (make-goldens.sh)
#   make update-goldens    the same, and write the goldens that are new or differ into ../data/real/
# All take MARKERS, SETS, CASES, FORMATS, THREADS, BRANCHES, GOLDEN_COMMIT, TRONKO_REAL_CACHE,
# TRONKO_REAL_MIRROR, TRONKO_REAL_WORK, e.g. make test-real MARKERS="mifish fwh" SETS="dev5k dev heldout"
test-real:
	bash test-real.sh

verify-branches:
	bash verify-branches.sh

goldens:
	bash make-goldens.sh

update-goldens:
	GOLDENS_UPDATE=1 bash make-goldens.sh

.PHONY: test-real verify-branches goldens update-goldens
