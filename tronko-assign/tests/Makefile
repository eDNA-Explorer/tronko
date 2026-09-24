# Tests of tronko-assign (README.md). The unit tests build independently of the main build.
#   make test       build and run every unit test with the production flags (gcc -O3
#                   -fno-stack-protector, no -march), then tronko-assign's own `make test`
#                   (run_tests.sh: the golden tests need tronko-assign built)
#   make test-san   the same under AddressSanitizer and UndefinedBehaviorSanitizer
#   make clean      remove the test binaries
CC = gcc
CFLAGS = -O3 -fno-stack-protector
SANFLAGS = -O1 -g -fsanitize=address,undefined -fno-sanitize-recover=undefined -fno-omit-frame-pointer
WFLAGS = -Wall -Wextra
KSW_WFLAGS = -Wall
BWA = ../bwa_source_files
KSW_DATA = data/ksw_fixture_calls.kswd.gz

test: test_vote_tally test_ksw_extend
	./test_vote_tally
	./test_ksw_extend $(KSW_DATA)
	$(MAKE) -C .. test

test-san: test_vote_tally_san test_ksw_extend_san
	ASAN_OPTIONS=detect_leaks=1 ./test_vote_tally_san
	ASAN_OPTIONS=detect_leaks=1 ./test_ksw_extend_san $(KSW_DATA)

test_vote_tally: unit/test_vote_tally.c ../vote_tally.h
	$(CC) $(CFLAGS) $(WFLAGS) -I.. -o $@ unit/test_vote_tally.c

test_vote_tally_san: unit/test_vote_tally.c ../vote_tally.h
	$(CC) $(SANFLAGS) $(WFLAGS) -I.. -o $@ unit/test_vote_tally.c

test_ksw_extend: unit/test_ksw_extend.c $(BWA)/ksw.c $(BWA)/ksw.h
	$(CC) $(CFLAGS) $(KSW_WFLAGS) -I$(BWA) -o $@ unit/test_ksw_extend.c $(BWA)/ksw.c -lz -lm

test_ksw_extend_san: unit/test_ksw_extend.c $(BWA)/ksw.c $(BWA)/ksw.h
	$(CC) $(SANFLAGS) $(KSW_WFLAGS) -I$(BWA) -o $@ unit/test_ksw_extend.c $(BWA)/ksw.c -lz -lm

clean:
	$(RM) test_vote_tally test_vote_tally_san test_ksw_extend test_ksw_extend_san

.PHONY: test test-san clean
