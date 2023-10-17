// Author: Martin C. Frith 2020
// SPDX-License-Identifier: GPL-3.0-or-later

// Code to scan a DNA sequence, looking for certain "words".

// First, call set().  Give it numOfPatterns "patterns", each of
// length wordLength, concatenated and pointed to by dnaMatches.  Each
// position in a pattern should indicate a subset of {a=1, c=2, g=4,
// t=8}, by a sum of {1,2,4,8}.

// numOfPatterns must be <= dnaWordsFinderNull.  If any 2 patterns can
// match the same DNA sequence, set() throws an error.

// wordLength must be small enough that the number of possible words
// fits in an unsigned int.

// We assume that a DNA sequence is encoded as a string of numbers:
// lettersToNumbers should map ASCII ACGTacgt to those numbers.

// You can then use count(), to count the number of occurrences of
// each pattern in a sequence.

// Or you can scan along a sequence, detecting each pattern, by using
// init() and next(), just as is done by count().

// next() returns dnaWordsFinderNull at positions that don't match any
// pattern, else it returns a number between 0 and numOfPatterns-1,
// indicating which pattern.

#ifndef DNA_WORDS_FINDER_HH
#define DNA_WORDS_FINDER_HH

#include <stddef.h>

#include <vector>

const unsigned dnaWordsFinderNull = 255;

struct DnaWordsFinder {
  typedef unsigned char uchar;

  unsigned wordLength;
  unsigned bitsPerBase;  // 1 for ry alphabet, 2 for acgt alphabet
  unsigned mask;
  unsigned numOfMatchedWords;
  std::vector<uchar> wordLookup;
  uchar baseToCode[256];

  void set(unsigned wordLength, unsigned numOfPatterns,
	   const uchar *dnaMatches, const uchar *lettersToNumbers,
	   bool isMaskLowercase);

  const uchar *init(const uchar *seqBeg, const uchar *seqEnd,
		    unsigned *hash) const {
    unsigned j = wordLength;
    while (j > 1 && seqBeg < seqEnd) {
      unsigned c = baseToCode[*seqBeg];
      ++seqBeg;
      if (c != dnaWordsFinderNull) {
	*hash = ((*hash) << bitsPerBase) + c;
	--j;
      } else {  // non-ACGT letter
	j = wordLength;
      }
    }
    return seqBeg;
  }

  unsigned next(unsigned *hash, unsigned newCode) const {
    *hash = ((*hash) << bitsPerBase) + newCode;
    return wordLookup[(*hash) & mask];
  }

  // Add counts of each word in [seqBeg, seqEnd) into "counts".
  void count(const uchar *seqBeg, const uchar *seqEnd, size_t *counts) const {
    unsigned hash = 0;
    seqBeg = init(seqBeg, seqEnd, &hash);
    while (seqBeg < seqEnd) {
      unsigned c = baseToCode[*seqBeg];
      ++seqBeg;
      if (c != dnaWordsFinderNull) {
	unsigned w = next(&hash, c);
	// xxx faster to count non-null words only?
	++counts[w];
      } else {  // non-ACGT letter
	seqBeg = init(seqBeg, seqEnd, &hash);
      }
    }
  }

private:
  bool isMatch(unsigned wordCode, const uchar *dnaMatches) {
    for (unsigned k = 0; k < wordLength; ++k) {
      unsigned baseCode = wordCode % (1 << bitsPerBase);
      unsigned matchCode = dnaMatches[wordLength - k - 1];
      if (!(matchCode & (1 << baseCode))) return false;
      wordCode >>= bitsPerBase;
    }
    return true;
  }
};

#endif
