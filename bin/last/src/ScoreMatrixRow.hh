// Copyright 2010, 2011 Martin C. Frith

// These definitions are for manipulating a score matrix, or a PSSM,
// as a 2-dimensional array.  The hope is that 2D arrays are extremely
// fast.  But arrays are "evil", so this may not be the best way...

#ifndef SCORE_MATRIX_ROW_HH
#define SCORE_MATRIX_ROW_HH

#include <limits.h>

namespace cbrc{

enum { scoreMatrixRowSize = ALPHABET_CAPACITY };

typedef int ScoreMatrixRow[scoreMatrixRowSize];

// Substitution score for delimiter symbols at the ends of sequences.
// It should be highly negative, to terminate alignments immediately,
// but not so negative that it causes overflow errors.

// The delimiter score when using short ints:
const int shortDelimiterScore = SHRT_MIN/2 + SCHAR_MIN;

// We want: short(delimiterScore) = shortDelimiterScore
const int delimiterScore = INT_MIN/2 + (unsigned short)shortDelimiterScore;

enum { INF = -delimiterScore };

}

#endif
