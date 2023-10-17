// Copyright 2008, 2009, 2011, 2014 Martin C. Frith

// This struct holds a score matrix for aligning pairs of residues,
// e.g. blosum62.  The delimiter symbol (space) aligned to anything
// gets a score of -INF.

// Maybe split this struct into two: ScoreMatrixEasy and ScoreMatrixFast?

#ifndef SCOREMATRIX_HH
#define SCOREMATRIX_HH

#include "ScoreMatrixRow.hh"

#include <stddef.h>

#include <iosfwd>
#include <string>
#include <vector>

namespace cbrc{

typedef unsigned char uchar;

struct ScoreMatrix{
  static const char *canonicalName( const std::string& name );
  static std::string stringFromName( const std::string& name );

  void setMatchMismatch(int matchScore,  // usually > 0
			int mismatchCost,  // usually > 0
			const std::string& symbols);  // case is preserved

  void fromString( const std::string& s );

  size_t numOfRows() const { return cells.size(); }
  size_t numOfCols() const { return cells[0].size(); }

  bool isCodonRows() const { return rowSymbols.size() > numOfRows(); }
  bool isCodonCols() const { return colSymbols.size() > numOfCols(); }

  void init(const uchar symbolToIndex[]);  // unspecified letters get minScore

  // Add scores for e.g. "W" meaning A or T.  The 2nd and 3rd
  // arguments specify whether to add a row/column for the
  // fully-ambiguous letter (N if DNA, else X).
  void addAmbiguousScores(bool isDna,
			  bool isFullyAmbiguousRow, bool isFullyAmbiguousCol,
			  const uchar symbolToIndex[],
			  double scale,  // "lambda" for getting probabilities
			  const double rowSymbolProbs[],
			  const double colSymbolProbs[]);

  void writeCommented( std::ostream& stream ) const;  // write preceded by "#"

  bool hasLetterFrequencies() const
  { return rowFrequencies.size() && colFrequencies.size(); }

  // store normalized letter frequencies in rowProbs and colProbs
  void calcLetterProbs(double *rowProbs, unsigned rowSize,
		       double *colProbs, unsigned colSize,
		       const uchar symbolToIndex[]) const;

  std::string rowSymbols;  // row headings (letters)
  std::string colSymbols;  // column headings (letters)
  std::vector< std::vector<int> > cells;  // scores
  std::vector<double> rowFrequencies;
  std::vector<double> colFrequencies;
  int caseSensitive[scoreMatrixRowSize][scoreMatrixRowSize];
  int caseInsensitive[scoreMatrixRowSize][scoreMatrixRowSize];
  int minScore;
  int maxScore;
};

std::istream& operator>>( std::istream& stream, ScoreMatrix& mat );
std::ostream& operator<<( std::ostream& stream, const ScoreMatrix& mat );

}

#endif
