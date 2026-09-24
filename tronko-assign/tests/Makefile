# Unit tests for tronko-assign, independent of the main build.
#   make test       build and run every unit test with the production flags (gcc -O3
#                   -fno-stack-protector, no -march), then the node-store unit test and the
#                   integration tests through tronko-assign's own `make test` (run_tests.sh;
#                   they need tronko-assign built, and report SKIP without their fixtures)
#   make test-san   the same under AddressSanitizer and UndefinedBehaviorSanitizer
#   make test-shadow  build tronko-assign with -DNW_SHADOW (the whole program, a minute or two)
#                   and run the repository fixtures with it (integration/test_nw_shadow.sh;
#                   TRONKO_FIXTURE_DIR as for test_assignment_production_parity.sh)
#   make clean      remove the test binaries
CC = gcc
CFLAGS = -O3 -fno-stack-protector
SANFLAGS = -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer
WFLAGS = -Wall -Wextra

NW_WFLAGS = -Wall
SRC = ..
NW_DEPS = $(SRC)/alignment.c $(SRC)/alignment.h $(SRC)/alignment_scoring.c \
          $(SRC)/alignment_scoring.h $(SRC)/alignment_macros.h $(SRC)/needleman_wunsch.c \
          $(SRC)/needleman_wunsch.h
NW_DATA = data/nw_fixture_alignments.txt.gz
# random alignments of test_nw_fill: the sanitizer build is about ten times slower
NW_RANDOM = 200000
NW_RANDOM_SAN = 20000

test: test_vote_tally test_nw_fill
	./test_vote_tally
	./test_nw_fill $(NW_DATA) $(NW_RANDOM)
	$(MAKE) -C .. test

test-san: test_vote_tally_san test_nw_fill_san
	ASAN_OPTIONS=detect_leaks=1 ./test_vote_tally_san
	ASAN_OPTIONS=detect_leaks=1 ./test_nw_fill_san $(NW_DATA) $(NW_RANDOM_SAN)

test-shadow: tronko-assign-nw-shadow
	TRONKO_ASSIGN_BIN=$(CURDIR)/tronko-assign-nw-shadow bash integration/test_nw_shadow.sh

test_vote_tally: unit/test_vote_tally.c ../vote_tally.h
	$(CC) $(CFLAGS) $(WFLAGS) -I.. -o $@ unit/test_vote_tally.c

test_vote_tally_san: unit/test_vote_tally.c ../vote_tally.h
	$(CC) $(SANFLAGS) $(WFLAGS) -I.. -o $@ unit/test_vote_tally.c

# unit/test_nw_fill.c includes alignment.c (its fills are static) and links the other two files
test_nw_fill: unit/test_nw_fill.c $(NW_DEPS)
	$(CC) $(CFLAGS) $(NW_WFLAGS) -I$(SRC) -o $@ unit/test_nw_fill.c $(SRC)/alignment_scoring.c $(SRC)/needleman_wunsch.c -lz

# alignment_scoring.c, unchanged seq-align code, gets AddressSanitizer only: where char is signed
# (x86-64) its scoring_lookup() indexes arrays with negative values for bytes above 0x7F, which
# UndefinedBehaviorSanitizer reports, and the test feeds it such bytes on purpose (see the test)
test_nw_fill_san: unit/test_nw_fill.c $(NW_DEPS)
	$(CC) $(SANFLAGS) -fno-sanitize=undefined $(NW_WFLAGS) -I$(SRC) -c -o alignment_scoring_san.o $(SRC)/alignment_scoring.c
	$(CC) $(SANFLAGS) $(NW_WFLAGS) -I$(SRC) -o $@ unit/test_nw_fill.c alignment_scoring_san.o $(SRC)/needleman_wunsch.c -lz
	$(RM) alignment_scoring_san.o

# the main build's compiler line (tronko-assign/Makefile) with -DNW_SHADOW; see nw_shadow.mk
tronko-assign-nw-shadow: $(NW_DEPS)
	$(MAKE) -C $(SRC) -f Makefile -f tests/nw_shadow.mk nw-shadow

clean:
	$(RM) test_vote_tally test_vote_tally_san unit/test_nodestore
	$(RM) test_nw_fill test_nw_fill_san tronko-assign-nw-shadow

.PHONY: test test-san test-shadow clean
