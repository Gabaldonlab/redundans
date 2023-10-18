// Copyright 2017 Martin C. Frith

#ifndef LAST_HH
#define LAST_HH

#include "Alphabet.hh"
#include "CyclicSubsetSeed.hh"
#include "MultiSequence.hh"
#include "SequenceFormat.hh"
#include "SubsetSuffixArray.hh"
#include "dna_words_finder.hh"
#include "qualityScoreUtil.hh"

namespace cbrc {

typedef MultiSequence::indexT indexT;

const size_t posSize = sizeof(PosPart) * posParts;
const size_t offSize = sizeof(OffPart) * offParts;

const size_t posLimit = size_t(-1) >> ((sizeof(size_t) - posSize) * CHAR_BIT);

inline void err(const char *s) { throw std::runtime_error(s); }

inline void throwSeqTooBig() { err("encountered a sequence that's too long"); }

inline void initSequences(MultiSequence &m, const Alphabet &a,
			  bool isTranslated, bool isAppendStopSymbol) {
  m.initForAppending(isTranslated ? 3 : 1, isAppendStopSymbol);
  a.tr(m.seqWriter(), m.seqWriter() + m.seqBeg(0));
}

inline void encodeSequences(MultiSequence &m, sequenceFormat::Enum f,
			    const Alphabet &a, bool isKeepLowercase,
			    indexT start) {
  size_t beg = m.seqBeg(start);
  size_t end = m.seqBeg(m.finishedSequences());
  a.tr(m.seqWriter() + beg, m.seqWriter() + end, isKeepLowercase);
  if (isPhred(f)) {
    checkQualityCodes(m.qualityReader() + beg,
		      m.qualityReader() + end, qualityOffset(f));
  }  // assumes one quality code per letter
}

// Read the next sequence, adding it to the MultiSequence
inline std::istream &appendSequence(MultiSequence &m, std::istream &in,
				    indexT maxSeqLen, sequenceFormat::Enum f,
				    const Alphabet &a, bool isMaskLowercase) {
  if (f == sequenceFormat::fasta) {
    m.appendFromFasta(in, maxSeqLen);
  } else if (f == sequenceFormat::fastx) {
    m.appendFromFastx(in, maxSeqLen, false);
  } else if (f == sequenceFormat::fastxKeep) {
    m.appendFromFastx(in, maxSeqLen, true);
  } else if (f == sequenceFormat::prb) {
    m.appendFromPrb(in, maxSeqLen, a.size, a.decode);
  } else if (f == sequenceFormat::pssm) {
    m.appendFromPssm(in, maxSeqLen, a.encode, isMaskLowercase);
  } else {
    m.appendFromFastq(in, maxSeqLen, true);
  }

  return in;
}

inline void makeWordsFinder(DnaWordsFinder &wordsFinder,
			    const CyclicSubsetSeed *seeds, size_t numOfSeeds,
			    const uchar *lettersToNumbers,
			    bool isMaskLowercase) {
  wordsFinder.wordLength = 0;
  size_t wordLength = maxRestrictedSpan(seeds, numOfSeeds);
  if (wordLength) {
    if (numOfSeeds > dnaWordsFinderNull) {
      err("I can't handle so many word-restricted seed patterns, sorry");
    }
    if (wordLength >= CHAR_BIT * sizeof(unsigned) / 2) {
      err("I can't handle such long word restrictions in seed patterns, sorry");
    }
    assert(numOfSeeds > 0);
    size_t lengthOfAllWords = numOfSeeds * wordLength;
    std::vector<uchar> dnaMatches(lengthOfAllWords);
    for (size_t i = 0; i < numOfSeeds; ++i) {
      seeds[i].matchingDna(&dnaMatches[i * wordLength], wordLength);
    }
    wordsFinder.set(wordLength, numOfSeeds, &dnaMatches[0],
		    lettersToNumbers, isMaskLowercase);
  }
}

}

#endif
