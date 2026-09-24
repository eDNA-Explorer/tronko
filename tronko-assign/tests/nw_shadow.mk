# Read after tronko-assign/Makefile, from tronko-assign/:
#   make -f Makefile -f tests/nw_shadow.mk nw-shadow
# Builds the plain tronko-assign with -DNW_SHADOW as tests/tronko-assign-nw-shadow, leaving
# tronko-assign/tronko-assign alone. In that build every alignment the two-pass fill computes
# is recomputed with the original fill, the three matrices are compared, and any difference
# aborts; at exit it prints "NW_SHADOW alignments_compared=<n> mismatches=0" to stderr.
nw-shadow:
	$(CC) $(OPTIMIZATION) $(ARCH_FLAGS) -DNW_SHADOW $(MATH_FLAGS) $(MEMOPT_FLAGS) $(STACKPROTECT) -o tests/tronko-assign-nw-shadow $(NEEDLEMANWUNSCH) $(HASHMAP) $(BWA) $(WFA2) $(SOURCES) $(LIBS) $(OPENMP_FLAGS)

.PHONY: nw-shadow
