// Copyright 2016 Martin C. Frith

// These routines extend an alignment in a given direction (forward or
// reverse) from given start points in two sequences.

// The algorithm is adapted from Section 3 of "A greedy algorithm for
// aligning DNA sequences" by Z Zhang, S Schwartz, L Wagner, W Miller,
// J Comput Biol. 2000 7(1-2):203-214.

// The start points point at the first positions we'll try to align.

// To use: first call "align", which calculates the alignment but only
// returns its score.  To get the actual alignment, call
// "getNextChunk" to get each gapless chunk.

// The sequences had better end with delimiter characters.

// The "scorer" indicates which letter pairs are considered matches.
// The algorithm uses a match score (taken from scorer[0][0]) and a
// mismatch cost (taken from scorer[0][1]).

// The match score must be divisible by 2.

#ifndef GREEDY_XDROP_ALIGNER_HH
#define GREEDY_XDROP_ALIGNER_HH

#include "ScoreMatrixRow.hh"

#include <stddef.h>
#include <vector>

namespace cbrc {

typedef unsigned char uchar;

class GreedyXdropAligner {
public:
  int align(const uchar *seq1,  // start point in the 1st sequence
	    const uchar *seq2,  // start point in the 2nd sequence
	    bool isForward,  // forward or reverse extension?
	    const ScoreMatrixRow *scorer,  // the substitution score matrix
	    int maxScoreDrop,
	    uchar delimiter);

  // Call this repeatedly to get each gapless chunk of the alignment.
  // The chunks are returned in far-to-near order.  The chunk's end
  // coordinates in each sequence (relative to the start of extension)
  // and length are returned in the 3 out-parameters.  If there are no
  // more chunks, the 3 parameters are unchanged and "false" is
  // returned.
  bool getNextChunk(size_t &end1, size_t &end2, size_t &length);

private:
  std::vector<size_t> rowOrigins;
  std::vector<int> furthest;  // analogous to the R array in Zhang et al.
  std::vector<int> minScore0s;  // analogous to the T array in Zhang et al.

  // Our position during the trace-back:
  int bestDistance;
  int bestDiagonal;

  const int *sources(int distance, int diagonal) const
  { return &furthest[rowOrigins[distance] + diagonal]; }
};

}

#endif
