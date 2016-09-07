// Copyright 2015 Martin C. Frith

// This class masks simple repeats (e.g. acacacacacac) in sequences,
// in one of 3 modes: normal mode, protein mode, or AT-rich DNA mode.

// The masking is insensitive to uppercase/lowercase of the input
// sequences.

// "isATrichDna" takes precedence over "isProtein".

// The "alphabet" (e.g. ACGT) is used only in normal mode.

// It is assumed that the sequences will be supplied after mapping
// letters to small integers (e.g. ACGT -> 0123).  "letterToIndex"
// specifies this mapping.  It is assumed that the small integers are
// less than scoreMatrixRowSize.

// "maskTable" defines how to do the masking: it maps small integers
// to masked integers.

#ifndef TANTAN_MASKER_HH
#define TANTAN_MASKER_HH

#include "ScoreMatrixRow.hh"
#include "tantan.hh"
#include <string>

namespace cbrc {

typedef unsigned char uchar;

class TantanMasker {
public:
  void init(bool isProtein,
	    bool isATrichDna,
	    const std::string &alphabet,
	    const uchar *letterToIndex);

  void mask(uchar *seqBeg, uchar *seqEnd, const uchar *maskTable) const {
    tantan::maskSequences(seqBeg, seqEnd, maxRepeatOffset, probMatrixPointers,
			  repeatProb, 0.05, 0.9, 0, 0, 0.5, maskTable );
  }

private:
  int maxRepeatOffset;
  double repeatProb;
  double probMatrix[scoreMatrixRowSize][scoreMatrixRowSize];
  double *probMatrixPointers[scoreMatrixRowSize];
};

}

#endif
