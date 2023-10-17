// Copyright 2009 Toshiyuki Sato

// Setup:
// First, read a genetic code with either operator>> or fromString.
// Then, initialize with either codeTableSet or initCodons.

#ifndef GENETICCODE_HH
#define GENETICCODE_HH

#include <string>
#include <vector>
#include <iosfwd>
#include <cassert>
#include <stddef.h>  // size_t

namespace cbrc{

typedef unsigned char uchar;

class Alphabet;

const int maxDnaAlphabetSize = 54;

class GeneticCode{
 private:
  std::string AAs;
  std::string Base[3];
  std::vector<uchar> genome2residue;
  std::vector<uchar> genome2residueWithoutLowercaseMasking;
  uchar codonToAminoAcid[256];

  static int codon2number(const uchar *codon) {
    const int n = maxDnaAlphabetSize;
    return n * n * codon[0] + n * codon[1] + codon[2];
  }

  static int codon2number2( const uchar* codon, const Alphabet& dnaAlph );

  friend std::istream& operator>>( std::istream& stream, GeneticCode& codon );

 public:
  // Converts a name to a text string defining a genetic code.
  // If the name isn't known, it assumes it's a file and tries to read it.
  static std::string stringFromName(const std::string &name);

  void fromString( const std::string& s );

  // Setup translation from DNA to amino acids
  void codeTableSet( const Alphabet& aaAlph, const Alphabet& dnaAlph );

  // Setup translation from DNA to codons (0=aaa, 1=aac, ..., 63=ttt,
  // 64=delimiter, 65=unknown).  Also setup codonToAminoAcid.  Any DNA
  // triplet with non-ACGT (or with lowercase if isMaskLowercase is
  // true) gets translated to 65=unknown.
  // If isMaskLowercase and isUnmaskLowercase are both true, also
  // setup translateWithoutMasking.
  void initCodons( const uchar *ntToNumber, const uchar *aaToNumber,
		   bool isMaskLowercase, bool isUnmaskLowercase );

  void translate( const uchar* beg, const uchar* end, uchar* dest ) const;

  void translateWithoutMasking( const uchar* beg, const uchar* end,
				uchar* dest ) const;

  uchar translation( const uchar* codon ) const
  { return genome2residue[ codon2number( codon ) ]; }

  const uchar *getCodonToAmino() const { return codonToAminoAcid; }
};

// Convert an amino-acid (translated) coordinate to a DNA coordinate
inline size_t aaToDna( size_t aaCoordinate, size_t frameSize ){
  if( frameSize == 0 ) return aaCoordinate;  // for non-translated sequences
  size_t frame = aaCoordinate / frameSize;
  size_t offset = aaCoordinate % frameSize;
  return frame + offset * 3;
}

// Convert a DNA coordinate to an amino-acid (translated) coordinate
inline size_t dnaToAa( size_t dnaCoordinate, size_t frameSize ){
  if( frameSize == 0 ) return dnaCoordinate;  // for non-translated sequences
  size_t frame = dnaCoordinate % 3;
  size_t offset = dnaCoordinate / 3;
  return frame * frameSize + offset;
}

// Convert begin and end coordinates to a size and a frameshift
inline void sizeAndFrameshift( size_t beg, size_t end,
			       size_t frameSize,  // 0 means not translated
			       size_t& size, size_t& frameshift ){
  if( frameSize ){  // if it's a translated sequence:
    size_t dnaBeg = aaToDna( beg, frameSize );
    size_t dnaEnd = aaToDna( end, frameSize );
    size_t dnaSize = dnaEnd - dnaBeg;
    assert( dnaBeg <= dnaEnd + 1 );  // allow a -1 frameshift
    size = ( dnaSize + 1 ) / 3;
    frameshift = ( dnaSize + 3 ) % 3;
  }
  else{  // if it's not a translated sequence:
    assert( beg <= end );
    size = end - beg;
    frameshift = 0;
  }
}

} // end namespace cbrc

#endif
