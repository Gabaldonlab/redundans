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

// If the isMaskLowercase argument of the reading routines is true,
// then all lowercase letters will get mapped to the DELIMITER subset,
// otherwise they will be treated like their uppercase equivalents.

#ifndef CYCLIC_SUBSET_SEED_HH
#define CYCLIC_SUBSET_SEED_HH

#include <iosfwd>
#include <string>
#include <vector>

namespace cbrc{

typedef unsigned char uchar;

class CyclicSubsetSeed{
public:
  enum { MAX_LETTERS = 64 };
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

  // Reads lines from "in" until it finds a pattern line.  Any seed
  // alphabet lines are appended to "seedAlphabet".  If it finds a
  // pattern line, it stores it in "pattern" and returns true, else it
  // returns false.  Blank lines and comment lines starting with # are
  // ignored.
  static bool nextPattern( std::istream& in,
			   std::vector< std::string >& seedAlphabet,
			   std::string& pattern );

  void clear() { subsetLists.clear(); subsetMaps.clear(); }

  // Seed letters are case-sensitive.
  void init( const std::vector< std::string >& seedAlphabet,
	     const std::string& pattern,
	     bool isMaskLowercase,
	     const uchar letterCode[] );

  // "inputLine" should be a grouping of sequence letters.
  void appendPosition( std::istream& inputLine,
		       bool isMaskLowercase,
		       const uchar letterCode[] );

  // Writes the grouping of sequence letters at the given position.
  // The position must be less than the span.
  void writePosition( std::ostream& out, unsigned position ) const;

  unsigned span() const{
    return subsetLists.size();
  }

  const uchar* subsetMap( unsigned depth ) const{
    return &subsetMaps[0] + (depth % span()) * MAX_LETTERS;
  }

  unsigned subsetCount( unsigned depth ) const{
    return subsetLists[ depth % span() ].size();
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

private:
  std::vector< std::vector<std::string> > subsetLists;
  std::vector<uchar> subsetMaps;
};

}

#endif
