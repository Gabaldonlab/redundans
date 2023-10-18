// Copyright 2008, 2009, 2010, 2011, 2014 Martin C. Frith

#include "ScoreMatrix.hh"
#include "ScoreMatrixData.hh"
#include "qualityScoreUtil.hh"
#include "zio.hh"

#include <assert.h>
#include <ctype.h>
#include <string.h>

#include <algorithm>  // min, max
#include <iomanip>
#include <sstream>
#include <stdexcept>

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

static void makeUppercase(std::string& s) {
  for (size_t i = 0; i < s.size(); ++i) {
    unsigned char c = s[i];
    s[i] = toupper(c);
  }
}

namespace cbrc{

typedef std::runtime_error Err;

static int baseToNumber(char uppercaseBase) {
  switch (uppercaseBase) {
  case 'A': return 0;
  case 'C': return 1;
  case 'G': return 2;
  case 'T': return 3;
  case 'N': return 4;
  default: throw Err("bad codon in score matrix");
  }
}

static int codonToNumber(const char *uppercaseCodon) {
  int x = baseToNumber(uppercaseCodon[0]);
  int y = baseToNumber(uppercaseCodon[1]);
  int z = baseToNumber(uppercaseCodon[2]);
  if (x < 4 && y < 4 && z < 4) return x * 16 + y * 4 + z;
  if (x > 3 && y > 3 && z > 3) return 65;
  throw Err("bad codon in score matrix");
}

const char *ScoreMatrix::canonicalName( const std::string& name ){
  for( size_t i = 0; i < COUNTOF(scoreMatrixNicknames); ++i )
    if( name == scoreMatrixNicknames[i].nickname )
      return scoreMatrixNicknames[i].realname;
  return name.c_str();
}

std::string ScoreMatrix::stringFromName( const std::string& name ){
  std::string n = canonicalName( name );

  for( size_t i = 0; i < COUNTOF(scoreMatrices); ++i )
    if( n == scoreMatrices[i].name )
      return scoreMatrices[i].text;

  return slurp( n.c_str() );
}

void ScoreMatrix::setMatchMismatch(int matchScore, int mismatchCost,
				   const std::string& symbols) {
  rowSymbols.assign(symbols.begin(), symbols.end());
  colSymbols.assign(symbols.begin(), symbols.end());

  size_t size = symbols.size();
  cells.resize(size);

  for (size_t i = 0; i < size; ++i) {
    cells[i].assign(size, -mismatchCost);
    cells[i][i] = matchScore;
  }
}

void ScoreMatrix::fromString( const std::string& matString ){
  std::istringstream iss(matString);
  iss >> *this;
  if (!iss) throw Err("can't read the score matrix");
}

static unsigned s2i(const uchar symbolToIndex[], uchar c) {
  return symbolToIndex[c];
}

// index in slow matrix => index in fast matrix
static unsigned fastIndex(size_t slowIndex, const char symbols[],
			  bool isCodons, const uchar symbolToIndex[]) {
  return isCodons ? codonToNumber(symbols + slowIndex * 3)
    :               s2i(symbolToIndex, symbols[slowIndex]);
}

static unsigned fastIndex(const char *symbolPtr,
			  bool isCodons, const uchar symbolToIndex[]) {
  return isCodons ? codonToNumber(symbolPtr) : s2i(symbolToIndex, *symbolPtr);
}

static void upperAndLowerIndex(unsigned &upper, unsigned &lower,
			       size_t slowIndex, const char symbols[],
			       bool isCodons, const uchar symbolToIndex[],
			       unsigned tooBig) {
  if (isCodons) {
    upper = lower = codonToNumber(symbols + slowIndex * 3);
    assert(upper < tooBig);
  } else {
    char symbol = symbols[slowIndex];
    uchar s = symbol;
    upper = symbolToIndex[s];
    lower = symbolToIndex[tolower(s)];
    if (upper >= tooBig || lower >= tooBig) {
      throw Err(std::string("bad letter in score matrix: ") + symbol);
    }
  }
}

void ScoreMatrix::init(const uchar symbolToIndex[]) {
  unsigned fastMatrixSize = scoreMatrixRowSize;
  assert(!rowSymbols.empty());
  assert(!colSymbols.empty());

  makeUppercase(rowSymbols);
  makeUppercase(colSymbols);

  minScore = maxScore = cells[0][0];
  for (size_t i = 0; i < numOfRows(); ++i) {
    for (size_t j = 0; j < numOfCols(); ++j) {
      minScore = std::min( minScore, cells[i][j] );
      maxScore = std::max( maxScore, cells[i][j] );
    }
  }

  // set default score = minScore:
  for (unsigned i = 0; i < fastMatrixSize; ++i) {
    for (unsigned j = 0; j < fastMatrixSize; ++j) {
      caseSensitive[i][j] = minScore;
      caseInsensitive[i][j] = minScore;
    }
  }

  for (size_t i = 0; i < numOfRows(); ++i) {
    for (size_t j = 0; j < numOfCols(); ++j) {
      unsigned iu, il, ju, jl;
      upperAndLowerIndex(iu, il, i, rowSymbols.c_str(), isCodonRows(),
			 symbolToIndex, fastMatrixSize);
      upperAndLowerIndex(ju, jl, j, colSymbols.c_str(), isCodonCols(),
			 symbolToIndex, fastMatrixSize);
      caseSensitive[iu][jl] = std::min( cells[i][j], 0 );
      caseSensitive[il][ju] = std::min( cells[i][j], 0 );
      caseSensitive[il][jl] = std::min( cells[i][j], 0 );
      caseSensitive[iu][ju] = cells[i][j];  // careful: maybe il==iu or jl==ju
      caseInsensitive[iu][ju] = cells[i][j];
      caseInsensitive[iu][jl] = cells[i][j];
      caseInsensitive[il][ju] = cells[i][j];
      caseInsensitive[il][jl] = cells[i][j];
    }
  }

  // set a hugely negative score for the delimiter symbol:
  uchar delimiter = ' ';
  unsigned d = isCodonRows() ? 64 : symbolToIndex[delimiter];
  unsigned e = isCodonCols() ? 64 : symbolToIndex[delimiter];
  assert(d < fastMatrixSize);
  assert(e < fastMatrixSize);
  for (unsigned i = 0; i < fastMatrixSize; ++i) {
    caseSensitive[d][i] = -INF;
    caseSensitive[i][e] = -INF;
    caseInsensitive[d][i] = -INF;
    caseInsensitive[i][e] = -INF;
  }
}

static void calcSomeLetterProbs(double probs[], unsigned alphabetSizeForProbs,
				const std::vector<double> &freqs,
				const char symbols[], bool isCodons,
				const uchar symbolToIndex[]) {
  double sum = 0;
  for (size_t i = 0; i < freqs.size(); ++i) {
    unsigned j = fastIndex(i, symbols, isCodons, symbolToIndex);
    if (j < alphabetSizeForProbs) {
      if (freqs[i] < 0) throw Err("bad score matrix: letter frequency < 0");
      sum += freqs[i];
    }
  }
  if (sum <= 0) throw Err("bad score matrix: no positive letter frequencies");

  std::fill_n(probs, alphabetSizeForProbs, 0.0);

  for (size_t i = 0; i < freqs.size(); ++i) {
    unsigned j = fastIndex(i, symbols, isCodons, symbolToIndex);
    if (j < alphabetSizeForProbs) {
      probs[j] = freqs[i] / sum;
    }
  }
}

void ScoreMatrix::calcLetterProbs(double rowProbs[], unsigned rowSize,
				  double colProbs[], unsigned colSize,
				  const uchar symbolToIndex[]) const {
  calcSomeLetterProbs(rowProbs, rowSize, rowFrequencies, rowSymbols.c_str(),
		      isCodonRows(), symbolToIndex);
  calcSomeLetterProbs(colProbs, colSize, colFrequencies, colSymbols.c_str(),
		      isCodonCols(), symbolToIndex);
}

void ScoreMatrix::writeCommented( std::ostream& stream ) const{
  size_t symbolsPerRow = rowSymbols.size() / numOfRows();
  size_t symbolsPerCol = colSymbols.size() / numOfCols();
  size_t colWidth = (numOfCols() < 20) ? 3 : 2;
  if (colWidth < symbolsPerCol) colWidth = symbolsPerCol;

  stream << "# " << std::setw(symbolsPerRow) << "";
  for (size_t i = 0; i < numOfCols(); ++i) {
    stream << ' ' << std::setw(colWidth - symbolsPerCol) << "";
    stream.write(colSymbols.c_str() + i * symbolsPerCol, symbolsPerCol);
  }
  stream << '\n';

  for (size_t i = 0; i < numOfRows(); ++i) {
    stream << "# ";
    stream.write(rowSymbols.c_str() + i * symbolsPerRow, symbolsPerRow);
    for (size_t j = 0; j < numOfCols(); ++j) {
      stream << ' ' << std::setw(colWidth) << cells[i][j];
    }
    stream << '\n';
  }
}

std::istream& operator>>( std::istream& stream, ScoreMatrix& m ){
  std::string tmpRowSymbols;
  std::string tmpColSymbols;
  std::vector< std::vector<int> > tmpCells;
  std::vector<double> tmpRowFreqs;
  std::vector<double> tmpColFreqs;
  std::string line, word;
  size_t symbolsPerRow = 1;
  size_t symbolsPerCol = 1;

  while (stream) {
    if (!getline(stream, line)) {
      if (stream.eof() && !stream.bad()) stream.clear(std::ios::eofbit);
      break;
    }
    std::istringstream iss(line);
    if (!(iss >> word)) continue;  // skip blank lines
    if (tmpColSymbols.empty()) {
      if (word[0] == '#') continue;  // skip comment lines at the top
      if (word.size() == 3) symbolsPerCol = 3;
      do {
	if (word.size() != symbolsPerCol) stream.setstate(std::ios::failbit);
	tmpColSymbols.insert(tmpColSymbols.end(), word.begin(), word.end());
      } while (iss >> word);
    } else {
      size_t numOfColumns = tmpColSymbols.size() / symbolsPerCol;
      std::vector<int> row;
      int score;
      double freq;
      for (size_t i = 0; i < numOfColumns; ++i) {
	iss >> score;
	row.push_back(score);
      }
      if (tmpCells.empty() && word.size() == 3) symbolsPerRow = 3;
      if (word.size() == symbolsPerRow && iss) {
	tmpRowSymbols.insert(tmpRowSymbols.end(), word.begin(), word.end());
	tmpCells.push_back(row);
	if (iss >> freq) {
	  tmpRowFreqs.push_back(freq);
	  if (tmpRowFreqs.size() < tmpCells.size()) {
	    stream.setstate(std::ios::failbit);
	  }
	}
      } else {
	std::istringstream iss2(line);
	while (iss2 >> freq) {
	  tmpColFreqs.push_back(freq);
	}
	if (tmpColFreqs.size() > numOfColumns || tmpColFreqs.empty()) {
	  stream.setstate(std::ios::failbit);
	}
	break;
      }
    }
  }

  if (tmpCells.empty() || tmpRowFreqs.empty() != tmpColFreqs.empty()) {
    stream.setstate(std::ios::failbit);
  }

  if (stream) {
    m.rowSymbols.swap(tmpRowSymbols);
    m.colSymbols.swap(tmpColSymbols);
    m.cells.swap(tmpCells);
    m.rowFrequencies.swap(tmpRowFreqs);
    m.colFrequencies.swap(tmpColFreqs);
  }

  return stream;
}

const char *ntAmbiguities[] = {
  "M" "AC",
  "S" "CG",
  "K" "GT",
  "W" "TA",
  "R" "AG",
  "Y" "CT",
  "B" "CGT",
  "D" "AGT",
  "H" "ACT",
  "V" "ACG",
  "N" "ACGT"
};

const char *aaAmbiguities[] = {
  "X" "ACDEFGHIKLMNPQRSTVWY"
};

static bool isIn(const std::string& s, char x) {
  return find(s.begin(), s.end(), x) != s.end();
}

static const char *ambiguityList(const char *ambiguities[],
				 size_t numOfAmbiguousSymbols,
				 const char *symbolPtr, size_t symbolLen,
				 char scratch[]) {
  for (size_t i = 0; i < numOfAmbiguousSymbols; ++i) {
    if (ambiguities[i][0] == *symbolPtr) return ambiguities[i] + 1;
  }
  memcpy(scratch, symbolPtr, symbolLen);
  return scratch;
}

static double symbolProbSum(bool isOneSymbol, const uchar symbolToIndex[],
			    const char *symbols, const double probs[]) {
  if (isOneSymbol) return 1;
  double p = 0;
  for (size_t i = 0; symbols[i]; ++i) {
    unsigned x = s2i(symbolToIndex, symbols[i]);
    p += probs[x];
  }
  return p;
}

static int jointScore(const uchar symbolToIndex[], int **fastMatrix,
		      double scale,
		      const double rowSymbolProbs[],
		      const double colSymbolProbs[],
		      const char *rSymbols, const char *cSymbols,
		      size_t symbolsPerRow, size_t symbolsPerCol) {
  bool isOneRowSymbol = (rSymbols[symbolsPerRow] == 0);
  bool isOneColSymbol = (cSymbols[symbolsPerCol] == 0);

  double p = 0;
  for (const char *i = rSymbols; *i; i += symbolsPerRow) {
    for (const char *j = cSymbols; *j; j += symbolsPerCol) {
      unsigned x = fastIndex(i, symbolsPerRow > 1, symbolToIndex);
      unsigned y = fastIndex(j, symbolsPerCol > 1, symbolToIndex);
      double r = isOneRowSymbol ? 1 : rowSymbolProbs[x];
      double c = isOneColSymbol ? 1 : colSymbolProbs[y];
      p += r * c * probFromScore(scale, fastMatrix[x][y]);
    }
  }

  double rowProbSum = symbolProbSum(isOneRowSymbol, symbolToIndex, rSymbols,
				    rowSymbolProbs);
  double colProbSum = symbolProbSum(isOneColSymbol, symbolToIndex, cSymbols,
				    colSymbolProbs);

  return scoreFromProb(scale, p / (rowProbSum * colProbSum));
}

void ScoreMatrix::addAmbiguousScores(bool isDna, bool isFullyAmbiguousRow,
				     bool isFullyAmbiguousCol,
				     const uchar symbolToIndex[],
				     double scale,
				     const double rowSymbolProbs[],
				     const double colSymbolProbs[]) {
  int *fastMatrix[scoreMatrixRowSize];
  std::copy(caseInsensitive, caseInsensitive + scoreMatrixRowSize, fastMatrix);

  if ((isFullyAmbiguousRow && isCodonRows()) ||
      (isFullyAmbiguousCol && isCodonCols()))
    throw Err("codon ambiguity not implemented");

  size_t symbolsPerRow = rowSymbols.size() / numOfRows();
  size_t symbolsPerCol = colSymbols.size() / numOfCols();

  char scratch[4] = {0};

  const char **ambiguities = isDna ? ntAmbiguities : aaAmbiguities;
  size_t n = isDna ? COUNTOF(ntAmbiguities) : COUNTOF(aaAmbiguities);
  size_t numOfAmbiguousRows = isCodonRows() ? 0 : n - 1 + isFullyAmbiguousRow;
  size_t numOfAmbiguousCols = isCodonCols() ? 0 : n - 1 + isFullyAmbiguousCol;

  for (size_t k = 0; k < numOfAmbiguousCols; ++k) {
    char ambiguousSymbol = ambiguities[k][0];
    if (isIn(colSymbols, ambiguousSymbol)) continue;
    colSymbols.push_back(ambiguousSymbol);
    for (size_t i = 0; i < numOfRows(); ++i) {
      const char *symbolPtr = rowSymbols.c_str() + i * symbolsPerRow;
      const char *rSymbols = ambiguityList(ambiguities, numOfAmbiguousRows,
					   symbolPtr, symbolsPerRow, scratch);
      const char *cSymbols = ambiguities[k] + 1;
      int s = jointScore(symbolToIndex, fastMatrix, scale,
			 rowSymbolProbs, colSymbolProbs, rSymbols, cSymbols,
			 symbolsPerRow, symbolsPerCol);
      cells[i].push_back(s);
    }
  }

  for (size_t k = 0; k < numOfAmbiguousRows; ++k) {
    char ambiguousSymbol = ambiguities[k][0];
    if (isIn(rowSymbols, ambiguousSymbol)) continue;
    rowSymbols.push_back(ambiguousSymbol);
    cells.resize(cells.size() + 1);
    for (size_t j = 0; j < numOfCols(); ++j) {
      const char *symbolPtr = colSymbols.c_str() + j * symbolsPerCol;
      const char *rSymbols = ambiguities[k] + 1;
      const char *cSymbols = ambiguityList(ambiguities, numOfAmbiguousCols,
					   symbolPtr, symbolsPerCol, scratch);
      int s = jointScore(symbolToIndex, fastMatrix, scale,
			 rowSymbolProbs, colSymbolProbs, rSymbols, cSymbols,
			 symbolsPerRow, symbolsPerCol);
      cells.back().push_back(s);
    }
  }
}

}
