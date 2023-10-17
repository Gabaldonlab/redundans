// Copyright 2010 Martin C. Frith

// Functions that find gapless X-drop alignments between two sequences.

#ifndef GAPLESS_XDROP_HH
#define GAPLESS_XDROP_HH

#include "ScoreMatrixRow.hh"
#include <stddef.h>

namespace cbrc {

typedef unsigned char uchar;

// Gets the maximum score for any gapless alignment starting at (seq1,
// seq2) and extending forwards.  The score is not allowed to drop by
// more than maxScoreDrop.  The sequences had better end with
// sentinels that have score < -maxScoreDrop.
// The score might suffer overflow, for huge sequences and/or huge
// scores.  If the function detects this (not guaranteed), it throws
// an exception.
int forwardGaplessXdropScore(const uchar *seq1,
                             const uchar *seq2,
                             const ScoreMatrixRow *scorer,
                             int maxScoreDrop);

// As above, but extending backwards.
int reverseGaplessXdropScore(const uchar *seq1,
                             const uchar *seq2,
                             const ScoreMatrixRow *scorer,
                             int maxScoreDrop);

// Return the endpoint in seq1 of the shortest alignment starting at
// (seq1, seq2) and extending forwards, that has score equal to
// "score".  This score should be the one that was found by
// forwardGaplessXdropScore.
const uchar *forwardGaplessXdropEnd(const uchar *seq1,
                                    const uchar *seq2,
                                    const ScoreMatrixRow *scorer,
                                    int score);

// As above, but extending backwards.
const uchar *reverseGaplessXdropEnd(const uchar *seq1,
                                    const uchar *seq2,
                                    const ScoreMatrixRow *scorer,
                                    int score);

// Check whether the gapless alignment starting at (seq1, seq2) and
// ending at seq1end is "optimal".  Here, "optimal" means: the
// alignment has no prefix with score <= 0, no suffix with score <= 0,
// and no region with score < -maxScoreDrop.
bool isOptimalGaplessXdrop(const uchar *seq1,
                           const uchar *seq1end,
                           const uchar *seq2,
                           const ScoreMatrixRow *scorer,
                           int maxScoreDrop);

// Returns the score, and sets the reverse and forward extension
// lengths, for a gapless "overlap" alignment starting at (seq1,
// seq2).  "Overlap" means that the alignment must extend, in each
// direction, until it hits a score <= -INF (presumably from a
// sentinel indicating a sequence end).  If the alignment would have
// any region with score < -maxScoreDrop, -INF is returned and the
// extension lengths are not set.
int gaplessXdropOverlap(const uchar *seq1,
			const uchar *seq2,
			const ScoreMatrixRow *scorer,
			int maxScoreDrop,
			size_t &reverseLength,
			size_t &forwardLength);

// Calculate the score of the gapless alignment starting at (seq1,
// seq2) and ending at seq1end.
int gaplessAlignmentScore(const uchar *seq1,
                          const uchar *seq1end,
                          const uchar *seq2,
                          const ScoreMatrixRow *scorer);

}

#endif
