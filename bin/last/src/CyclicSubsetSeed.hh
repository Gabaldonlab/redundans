// Copyright 2009, 2010, 2013, 2014 Martin C. Frith

// A "subset seed" covers a range of sequence.  The size of this range
// is called its span.  At each position, it maps letters (encoded as
// numbers) to subsets of the letters (encoded as numbers).  The
// mapping may or may not be different at different positions.

// Subset seeds are described in:
// G Kucherov et al. J Bioinform Comput Biol. 2006 4(2):553-69.

// "Cyclic" means that the seed can cover arbitrary-size ranges, by
// cyclically repeating.

// At each seed position, the subsets are defined by a grouping of
// sequence letters, such as "AG CT".  Any other sequence letters
// (such as "N" in this case) will get mapped to the special DELIMITER
// subset.

// An "exact" position is one with no grouped letters, e.g. "A C G T".

// A "restricted" position is one that omits letters from the main
// sequence alphabet, e.g. "A G".

// At a restricted position, the omitted main-alphabet letters are not
// actually mapped to the DELIMITER subset.  If the position is exact,
// then the omitted letters are mapped to separate subsets, else they
// are mapped to the same subset.

// If the isMaskLowercase argument of the reading routines is true,
// then all lowercase letters will get mapped to the DELIMITER subset,
// otherwise they will be treated like their uppercase equivalents.

#ifndef CYCLIC_SUBSET_SEED_HH
#define CYCLIC_SUBSET_SEED_HH

#include <stddef.h>

#include <iosfwd>
#include <stdexcept>
#include <string>
#include <vector>

namespace cbrc{

typedef unsigned char uchar;

class CyclicSubsetSeed{
public:
  enum { MAX_LETTERS = ALPHABET_CAPACITY };
  enum { DELIMITER = 255 };

  // Converts a name to a text string defining one or more seeds.
  // If the name isn't known, it assumes it's a file and tries to read it.
  static std::string stringFromName( const std::string& name );

  // Converts patterns to a text string defining one or more seeds.
  // "patterns" should be something like: "1110TT,1001T1".  The "1"s
  // are must-match positions, the "0"s are don't care positions, and
  // "T" or "t" allows transitions but not transversions.  For
  // consistency with YASS/Iedera, you can also use "#" for match, "@"
  // for transition, and "_" or "-" for don't care.  You can have
  // multiple seed patterns separated by commas.  "sequenceLetters"
  // indicates the letters that can occur in "1" and "0" positions.
  static std::string stringFromPatterns( const std::string& patterns,
					 const std::string& sequenceLetters );

  // Converts patterns to a text string defining one or more seeds.
  // "patterns" should be something like: "RYynN@,RyR@nN".  Uppercase
  // letters are must-match positions; lowercase letters are
  // mismatch-tolerant positions.  The possible letters are NRYACGT:
  // they only allow bases that match according to IUPAC ambiguity
  // codes.  Finally, "@" allows any match or transition.
  static std::string stringFromDnaPatterns(std::string patterns);

  // Read subset-seed patterns from text and add them to patterns.
  // The text should start with lines defining a seed alphabet,
  // followed by lines with seed patterns.  Blank lines and comment
  // lines starting with # are ignored.  Seed letters are
  // case-sensitive.
  static void addPatterns(std::vector<CyclicSubsetSeed> &patterns,
			  const std::string &text,
			  bool isMaskLowercase,
			  const uchar letterCode[],
			  const std::string &mainSequenceAlphabet);

  void swap( CyclicSubsetSeed& x ){
    subsetLists.swap( x.subsetLists );
    subsetMaps.swap( x.subsetMaps );
    originalSubsetMaps.swap( x.originalSubsetMaps );
    numOfSubsetsPerPosition.swap( x.numOfSubsetsPerPosition );
  }

  // "inputLine" should be a grouping of sequence letters.
  void appendPosition( std::istream& inputLine,
		       bool isMaskLowercase,
		       const uchar letterCode[],
		       const std::string& mainSequenceAlphabet );

  // E.g. if the seed maps amino acids to subsets, and initialMap maps
  // codons to amino acids, then this gives us codons=>subsets:
  void compose(const uchar *initialMap) {
    size_t size = subsetMaps.size();
    for (size_t i = 0; i < size; i += MAX_LETTERS) {
      for (size_t j = 0; j < MAX_LETTERS; ++j) {
	subsetMaps[i + j] = originalSubsetMaps[i + initialMap[j]];
      }
    }
  }

  // Get the original map (e.g. amino acids => subsets) at the same
  // seed position as the composed map (e.g. codons => subsets):
  const uchar *originalSubsetMap(const uchar *composedMap) const {
    return &originalSubsetMaps[0] + (composedMap - &subsetMaps[0]);
  }

  // Writes the grouping of sequence letters at the given position.
  // The position must be less than the span.
  void writePosition(std::ostream& out, size_t position) const;

  size_t span() const {
    return subsetLists.size();
  }

  const uchar* subsetMap(size_t depth) const {
    return &subsetMaps[0] + (depth % span()) * MAX_LETTERS;
  }

  unsigned restrictedSubsetCount(size_t depth) const {
    return subsetLists[depth].size();
  }

  unsigned unrestrictedSubsetCount(size_t depth) const {
    return numOfSubsetsPerPosition[depth % numOfSubsetsPerPosition.size()];
  }

  // Number of positions up to & including the rightmost restricted position
  size_t restrictedSpan() const {
    for (size_t i = subsetLists.size(); i > 0; --i) {
      if (subsetLists[i-1].size() < numOfSubsetsPerPosition[i-1]) return i;
    }
    return 0;
  }

  const uchar* firstMap() const{
    return &subsetMaps[0];
  }

  const uchar* nextMap( const uchar* x ) const{
    const uchar* y = x + MAX_LETTERS;
    if( y == &subsetMaps[0] + subsetMaps.size() )
      y = &subsetMaps[0];
    return y;
  }

  const uchar* prevMap( const uchar* x ) const{
    if( x == &subsetMaps.front() )
      x = &subsetMaps.back() + 1;
    return x - MAX_LETTERS;
  }

  // Checks whether text1 is lexicographically less than text2, when
  // applying the cyclic subset seed pattern to both.  The "startMap"
  // argument enables us to start in the middle of the seed pattern.
  bool isLess(const uchar *text1, const uchar *text2,
	      const uchar *startMap) const {
    while (true) {
      uchar x = startMap[*text1];
      uchar y = startMap[*text2];
      if (x != y) return x < y;
      if (x == DELIMITER) return false;
      ++text1;
      ++text2;
      startMap = nextMap(startMap);
    }
  }

  // Which DNA sequences, of length wordLength (which must be <= the
  // span), can match the start of the seed.  The output is written in
  // dnaMatches, which should point to wordLength zeros.
  void matchingDna(uchar *dnaMatches, unsigned wordLength) const {
    for (unsigned i = 0; i < wordLength; ++i) {
      for (unsigned j = 0; j < subsetLists[i].size(); ++j) {
	for (unsigned k = 0; k < subsetLists[i][j].size(); ++k) {
	  switch (subsetLists[i][j][k]) {
	  case 'A': dnaMatches[i] |= 1; break;
	  case 'C': dnaMatches[i] |= 2; break;
	  case 'G': dnaMatches[i] |= 4; break;
	  case 'T': dnaMatches[i] |= 8; break;
	  default:
	    throw std::runtime_error("I can't handle non-DNA in "
				     "word-restricted seeds, sorry");
	  }
	}
      }
    }
  }

private:
  std::vector< std::vector<std::string> > subsetLists;
  std::vector<uchar> subsetMaps;
  std::vector<uchar> originalSubsetMaps;
  std::vector<unsigned> numOfSubsetsPerPosition;
};

inline size_t maxRestrictedSpan(const CyclicSubsetSeed *seeds, size_t n) {
  size_t x = 0;
  for (size_t i = 0; i < n; ++i) {
    size_t y = seeds[i].restrictedSpan();
    if (y > x) x = y;
  }
  return x;
}

}

#endif
