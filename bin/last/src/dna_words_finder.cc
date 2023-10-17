// Author: Martin C. Frith 2020
// SPDX-License-Identifier: GPL-3.0-or-later

#include "dna_words_finder.hh"

#include <assert.h>
#include <string.h>

#include <stdexcept>

static unsigned getBitsPerBase(const unsigned char *dnaMatches, unsigned size){
  for (unsigned i = 0; i < size; ++i) {
    if (dnaMatches[i] % 5) return 2;  // can't use ry alphabet
  }
  return 1;
}

void DnaWordsFinder::set(unsigned wordLength, unsigned numOfPatterns,
			 const uchar *dnaMatches,
			 const uchar *lettersToNumbers,
			 bool isMaskLowercase) {
  this->wordLength = wordLength;
  if (wordLength == 0) return;
  assert(numOfPatterns <= dnaWordsFinderNull);
  unsigned totalPatternLength = wordLength * numOfPatterns;
  bitsPerBase = getBitsPerBase(dnaMatches, totalPatternLength);
  unsigned numOfAllWords = 1;
  numOfAllWords <<= (wordLength * bitsPerBase);
  mask = numOfAllWords - 1;
  wordLookup.assign(numOfAllWords, dnaWordsFinderNull);
  numOfMatchedWords = 0;

  for (unsigned i = 0; i < numOfAllWords; ++i) {
    for (unsigned j = 0; j < numOfPatterns; ++j) {
      if (isMatch(i, dnaMatches + j * wordLength)) {
	if (wordLookup[i] != dnaWordsFinderNull) {
	  using std::runtime_error;
	  throw runtime_error("2 word-restricted seed patterns match 1 word");
	}
	wordLookup[i] = j;
	++numOfMatchedWords;
      }
    }
  }

  memset(baseToCode, dnaWordsFinderNull, sizeof baseToCode);

  for (unsigned i = 0; i < 4; ++i) {
    unsigned c = i % (1 << bitsPerBase);
    baseToCode[lettersToNumbers["ACGT"[i]]] = c;
    if (!isMaskLowercase) baseToCode[lettersToNumbers["acgt"[i]]] = c;
  }
}
