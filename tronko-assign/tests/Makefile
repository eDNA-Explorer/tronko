# Unit tests for tronko-assign, independent of the main build.
#   make test       build and run every unit test with the production flags (gcc -O3
#                   -fno-stack-protector, no -march), then the node-store unit test and the
#                   integration tests through tronko-assign's own `make test` (run_tests.sh;
#                   they need tronko-assign built, and report SKIP without their fixtures)
#   make test-san   the same under AddressSanitizer and UndefinedBehaviorSanitizer
#   make test-tsan  the threaded ones under ThreadSanitizer (x86-64; if TSan stops with
#                   "unexpected memory mapping", run it as: setarch x86_64 -R make test-tsan)
#   make test-shadow  build tronko-assign with -DNW_SHADOW (the whole program, a minute or two)
#                   and run the repository fixtures with it (integration/test_nw_shadow.sh;
#                   TRONKO_FIXTURE_DIR as for test_assignment_production_parity.sh)
#   make clean      remove the test binaries
CC = gcc
CFLAGS = -O3 -fno-stack-protector
SANFLAGS = -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer
TSANFLAGS = -O1 -g -fsanitize=thread
WFLAGS = -Wall -Wextra
KSW_WFLAGS = -Wall
BWA = ../bwa_source_files
KSW_DATA = data/ksw_fixture_calls.kswd.gz
NW_WFLAGS = -Wall
SRC = ..
NW_DEPS = $(SRC)/alignment.c $(SRC)/alignment.h $(SRC)/alignment_scoring.c \
          $(SRC)/alignment_scoring.h $(SRC)/alignment_macros.h $(SRC)/needleman_wunsch.c \
          $(SRC)/needleman_wunsch.h
NW_DATA = data/nw_fixture_alignments.txt.gz
# random alignments of test_nw_fill: the sanitizer build is about ten times slower
NW_RANDOM = 200000
NW_RANDOM_SAN = 20000
DENSE_WFLAGS = -Wall
LIBS = -lz -lm -lpthread
FIXTURE_INDEX = ../../tronko-build/example_datasets/single_tree/Charadriiformes.fasta

# test_dense_sa: bwt.c under test (with warnings), and the BWA index-building code the test uses
# to make synthetic indexes (bwt_bwtupdate_core in bwtindex.c and what it links)
DENSE_SA_BWA = $(BWA)/bwtindex.c $(BWA)/utils.c $(BWA)/bntseq.c $(BWA)/bwt_gen.c $(BWA)/QSufSort.c \
	$(BWA)/is.c $(BWA)/rope.c $(BWA)/rle.c $(BWA)/kstring.c
DENSE_SA_DEPS = unit/test_dense_sa.c $(BWA)/bwt.c $(BWA)/bwt.h
CHAIN_WFLAGS = -Wall
# BWA sources a test links besides the one it includes; fastmap.c needs the rest of tronko-assign.
BWA_LINK = $(BWA)/bntseq.c $(BWA)/bwa.c $(BWA)/bwamem_pair.c $(BWA)/bwamem_extra.c $(BWA)/bwt.c \
	$(BWA)/ksw.c $(BWA)/kstring.c $(BWA)/kthread.c $(BWA)/utils.c
CHAIN_DATA = data/chain_flt_fixture.chf.gz
CHAIN_DEPS = $(wildcard $(BWA)/*.c $(BWA)/*.h) ../global.h

test: test_vote_tally test_ksw_extend test_nw_fill test_dense_sa test_chain_flt test_pos2rid
	./test_vote_tally
	./test_ksw_extend $(KSW_DATA)
	./test_nw_fill $(NW_DATA) $(NW_RANDOM)
	./test_dense_sa $(FIXTURE_INDEX)
	./test_chain_flt $(CHAIN_DATA)
	./test_pos2rid $(FIXTURE_INDEX)
	$(MAKE) -C .. test

test-san: test_vote_tally_san test_ksw_extend_san test_nw_fill_san test_dense_sa_san test_chain_flt_san test_pos2rid_san
	ASAN_OPTIONS=detect_leaks=1 ./test_vote_tally_san
	ASAN_OPTIONS=detect_leaks=1 ./test_ksw_extend_san $(KSW_DATA)
	ASAN_OPTIONS=detect_leaks=1 ./test_nw_fill_san $(NW_DATA) $(NW_RANDOM_SAN)
	ASAN_OPTIONS=detect_leaks=1 ./test_dense_sa_san $(FIXTURE_INDEX)
	ASAN_OPTIONS=detect_leaks=1 ./test_chain_flt_san $(CHAIN_DATA)
	ASAN_OPTIONS=detect_leaks=1 ./test_pos2rid_san $(FIXTURE_INDEX)

test-tsan: test_dense_sa_tsan
	TSAN_OPTIONS=halt_on_error=1 ./test_dense_sa_tsan $(FIXTURE_INDEX)

test-shadow: tronko-assign-nw-shadow
	TRONKO_ASSIGN_BIN=$(CURDIR)/tronko-assign-nw-shadow bash integration/test_nw_shadow.sh

test_vote_tally: unit/test_vote_tally.c ../vote_tally.h
	$(CC) $(CFLAGS) $(WFLAGS) -I.. -o $@ unit/test_vote_tally.c

test_vote_tally_san: unit/test_vote_tally.c ../vote_tally.h
	$(CC) $(SANFLAGS) $(WFLAGS) -I.. -o $@ unit/test_vote_tally.c

test_ksw_extend: unit/test_ksw_extend.c $(BWA)/ksw.c $(BWA)/ksw.h
	$(CC) $(CFLAGS) $(KSW_WFLAGS) -I$(BWA) -o $@ unit/test_ksw_extend.c $(BWA)/ksw.c -lz -lm

test_ksw_extend_san: unit/test_ksw_extend.c $(BWA)/ksw.c $(BWA)/ksw.h
	$(CC) $(SANFLAGS) $(KSW_WFLAGS) -I$(BWA) -o $@ unit/test_ksw_extend.c $(BWA)/ksw.c -lz -lm

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

# build_dense_sa <flags> <output>: the test and bwt.c with warnings, the rest of BWA without
define build_dense_sa
	$(CC) $(1) $(DENSE_WFLAGS) -I$(BWA) -c -o $(2).test.o unit/test_dense_sa.c
	$(CC) $(1) $(DENSE_WFLAGS) -I$(BWA) -c -o $(2).bwt.o $(BWA)/bwt.c
	$(CC) $(1) -w -o $(2) $(2).test.o $(2).bwt.o $(DENSE_SA_BWA) $(LIBS)
	$(RM) $(2).test.o $(2).bwt.o
endef

test_dense_sa: $(DENSE_SA_DEPS)
	$(call build_dense_sa,$(CFLAGS),$@)

test_dense_sa_san: $(DENSE_SA_DEPS)
	$(call build_dense_sa,$(SANFLAGS),$@)

test_dense_sa_tsan: $(DENSE_SA_DEPS)
	$(call build_dense_sa,$(TSANFLAGS),$@)

test_chain_flt: unit/test_chain_flt.c $(CHAIN_DEPS)
	$(CC) $(CFLAGS) $(CHAIN_WFLAGS) -o $@ unit/test_chain_flt.c $(filter-out $(BWA)/bwamem.c,$(BWA_LINK)) $(LIBS)

test_chain_flt_san: unit/test_chain_flt.c $(CHAIN_DEPS)
	$(CC) $(SANFLAGS) $(CHAIN_WFLAGS) -o $@ unit/test_chain_flt.c $(filter-out $(BWA)/bwamem.c,$(BWA_LINK)) $(LIBS)

test_pos2rid: unit/test_pos2rid.c $(CHAIN_DEPS)
	$(CC) $(CFLAGS) $(CHAIN_WFLAGS) -o $@ unit/test_pos2rid.c $(BWA)/bwamem.c $(BWA_LINK) $(LIBS)

test_pos2rid_san: unit/test_pos2rid.c $(CHAIN_DEPS)
	$(CC) $(SANFLAGS) $(CHAIN_WFLAGS) -o $@ unit/test_pos2rid.c $(BWA)/bwamem.c $(BWA_LINK) $(LIBS)

clean:
	$(RM) test_vote_tally test_vote_tally_san test_ksw_extend test_ksw_extend_san unit/test_nodestore
	$(RM) test_nw_fill test_nw_fill_san tronko-assign-nw-shadow
	$(RM) test_dense_sa test_dense_sa_san test_dense_sa_tsan
	$(RM) test_chain_flt test_chain_flt_san test_pos2rid test_pos2rid_san

.PHONY: test test-san test-tsan test-shadow clean
