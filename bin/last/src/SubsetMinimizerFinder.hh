// Copyright 2016 Martin C. Frith

// This class finds suffixes of a text that are "subset minimizers".

// Here, a minimizer is a suffix that is the lexicographic minimum of
// the suffixes that start in some range [max(i - window + 1, beg), i]
// for any i in [beg, end).

// "Subset" means that the lexicographic comparison is done on
// suffixes transformed according to a CyclicSubsetSeed.

// Usage:

// First call "init" to prepare for range [beg, end) in a given text.

// Then repeatedly call isMinimizer, to check whether the suffix
// starting at "pos" is a minimizer.  The suffixes must be checked in
// order: specifically, each value of pos must be >= beg and all
// previous values of pos.

#ifndef SUBSET_MINIMIZER_FINDER_HH
#define SUBSET_MINIMIZER_FINDER_HH

#include <vector>
#include <stddef.h>  // size_t

namespace cbrc {

typedef unsigned char uchar;

class CyclicSubsetSeed;

class SubsetMinimizerFinder {
public:
  void init(const CyclicSubsetSeed &seed,
	    const uchar *beg,
	    const uchar *end);

   bool isMinimizer(const CyclicSubsetSeed &seed,
		    const uchar *pos,
		    const uchar *end,
		    size_t window);

private:
  std::vector<const uchar *> minima;
};

}

#endif
