// Author: Martin C. Frith 2013
// SPDX-License-Identifier: GPL-3.0-or-later

#include "last-split.hh"
#include "mcf_last_splitter.hh"

#include <string.h>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <streambuf>

using namespace mcf;

class MyString {
public:
  MyString() : v(1), s(0), e(0) {}

  size_t size() const { return s; }

  char &operator[](size_t i) { return v[i]; }

  void resize(size_t i) {
    memmove(&v[0] + i, &v[0] + s, e - s);
    e -= s - i;
    s = i;
  }

  void erasePrefix(size_t len) {
    s -= len;
    e -= len;
    memmove(&v[0], &v[0] + len, e);
  }

  bool appendLine(std::istream &stream) {
    size_t i = s;
    for (;;) {
      char *x = static_cast<char *>(memchr(&v[0] + i, '\n', e - i));
      if (x) {
	*x = 0;
	s = x - &v[0] + 1;
	return true;
      }
      i = e;
      e += 256;  // xxx ???
      if (v.size() < e) v.resize(e);
      e = i + stream.rdbuf()->sgetn(&v[i], e - i);
      if (e == i) {
	if (i == s) return false;
	v[i] = 0;
	e = s = i + 1;
	return true;
      }
    }
  }

private:
  std::vector<char> v;
  size_t s;
  size_t e;
};

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static std::istream& openIn(const std::string& fileName, std::ifstream& ifs) {
  if (fileName == "-") return std::cin;
  ifs.open(fileName.c_str());
  if (!ifs) err("can't open file: " + fileName);
  return ifs;
}

// Does the string start with the prefix?
static bool startsWith(const char *s, const char *prefix) {
  for (;;) {
    if (*prefix == 0) return true;
    if (*prefix != *s) return false;
    ++s;
    ++prefix;
  }
}

// Does the string have no non-space characters?
static bool isBlankLine(const char *s) {
  for (;;) {
    if (*s == 0) return true;
    if (!std::isspace(*s)) return false;
    ++s;
  }
}

static bool isSpace(char c) {
  return c > 0 && c <= ' ';
}

static bool isSameName(const char *sLine1, const char *sLine2) {
  do { ++sLine1; } while (isSpace(*sLine1));
  do { ++sLine2; } while (isSpace(*sLine2));
  for (;;) {
    if (*sLine1 > ' ') {
      if (*sLine2 != *sLine1) return false;
    } else {
      return *sLine2 <= ' ';
    }
    ++sLine1;
    ++sLine2;
  }
}

static void transpose(std::vector< std::vector<int> > &matrix) {
  size_t rows = matrix.size();
  size_t cols = matrix[0].size();
  std::vector< std::vector<int> > m(cols);
  for (size_t i = 0; i < cols; ++i) {
    m[i].resize(rows);
    for (size_t j = 0; j < rows; ++j) {
      m[i][j] = matrix[j][i];
    }
  }
  m.swap(matrix);
}

static void doOneBatch(MyString &inputText,
		       const std::vector<size_t> &lineEnds,
		       const std::vector<unsigned> &mafEnds,
		       LastSplitter &splitter, const LastSplitOptions &opts,
		       const cbrc::SplitAlignerParams &params,
		       bool isAlreadySplit) {
  std::vector<char *> linePtrs(lineEnds.size());
  for (size_t i = 0; i < lineEnds.size(); ++i) {
    linePtrs[i] = &inputText[0] + lineEnds[i];
  }

  splitter.reserve(mafEnds.size() - 1);  // saves memory: no excess capacity
  for (unsigned i = 1; i < mafEnds.size(); ++i) {
    splitter.addMaf(&linePtrs[0] + mafEnds[i-1],
		    &linePtrs[0] + mafEnds[i], opts.isTopSeqQuery);
  }

  splitter.split(opts, params, isAlreadySplit);

  if (!splitter.isOutputEmpty()) {
    splitter.printOutput();
    splitter.clearOutput();
  }
}

static void addMaf(std::vector<unsigned> &mafEnds,
		   const std::vector<size_t> &lineEnds) {
  if (lineEnds.size() - 1 > mafEnds.back())  // if we have new maf lines:
    mafEnds.push_back(lineEnds.size() - 1);
}

static void eraseOldInput(MyString &inputText,
			  std::vector<size_t> &lineEnds,
			  std::vector<unsigned> &mafEnds) {
  size_t numOfOldLines = mafEnds.back();
  size_t numOfOldChars = lineEnds[numOfOldLines];
  inputText.erasePrefix(numOfOldChars);
  lineEnds.erase(lineEnds.begin(), lineEnds.begin() + numOfOldLines);
  for (size_t i = 0; i < lineEnds.size(); ++i) {
    lineEnds[i] -= numOfOldChars;
  }
  mafEnds.resize(1);
}

void lastSplit(LastSplitOptions& opts) {
  cbrc::SplitAlignerParams params;
  LastSplitter splitter;
  std::vector< std::vector<int> > scoreMatrix;
  std::string rowNames, colNames;
  std::string word, name, key;
  int state = 0;
  int sequenceFormat = 1;  // xxx ???
  int gapExistenceCost = -1;
  int gapExtensionCost = -1;
  int insExistenceCost = -1;
  int insExtensionCost = -1;
  int lastalScoreThreshold = -1;
  double scale = 0;
  double genomeSize = 0;
  MyString inputText;
  std::vector<size_t> lineEnds(1);  // offsets in inputText of line starts/ends
  std::vector<unsigned> mafEnds(1);  // which lines are in which MAF block
  unsigned sLineCount = 0;
  size_t qNameLineBeg = 0;
  bool isAlreadySplit = false;  // has the input already undergone last-split?

  for (unsigned i = 0; i < opts.inputFileNames.size(); ++i) {
    std::ifstream inFileStream;
    std::istream& input = openIn(opts.inputFileNames[i], inFileStream);
    while (inputText.appendLine(input)) {
      const char *linePtr = &inputText[0] + lineEnds.back();
      if (state == -1) {  // we are reading the score matrix within the header
	std::istringstream ls(linePtr);
	std::vector<int> row;
	int score;
	ls >> word >> name;
	while (ls >> score) row.push_back(score);
	if (word == "#" && name.size() == 1 &&
	    row.size() == colNames.size() && ls.eof()) {
	  rowNames.push_back(std::toupper(name[0]));
	  scoreMatrix.push_back(row);
	} else {
	  state = 0;
	}
      }
      if (state == 0) {  // we are reading the header
	std::istringstream ls(linePtr);
	std::string names;
	ls >> word;
	while (ls >> name) {
	  if (name.size() == 1) names.push_back(std::toupper(name[0]));
	  else break;
	}
	if (word == "#" && !names.empty() && !ls && scoreMatrix.empty()) {
	  colNames = names;
	  state = -1;
	} else if (linePtr[0] == '#') {
	  std::istringstream ls(linePtr);
	  while (ls >> word) {
	    std::istringstream ws(word);
	    getline(ws, key, '=');
	    if (key == "a") ws >> gapExistenceCost;
	    if (key == "b") ws >> gapExtensionCost;
	    if (key == "A") ws >> insExistenceCost;
	    if (key == "B") ws >> insExtensionCost;
	    if (key == "e") ws >> lastalScoreThreshold;
	    if (key == "t") ws >> scale;
	    if (key == "Q") ws >> sequenceFormat;
	    if (key == "letters") ws >> genomeSize;
	  }
	  // try to determine if last-split was already run (fragile):
	  if (startsWith(linePtr, "# m=")) isAlreadySplit = true;
	} else if (!isBlankLine(linePtr)) {
	  if (scoreMatrix.empty())
	    err("I need a header with score parameters");
	  if (gapExistenceCost < 0 || gapExtensionCost < 0 ||
	      insExistenceCost < 0 || insExtensionCost < 0 ||
	      lastalScoreThreshold < 0 || scale <= 0 || genomeSize <= 0)
	    err("can't read the header");
	  opts.setUnspecifiedValues(lastalScoreThreshold, scale);
	  if (opts.isTopSeqQuery) {
	    transpose(scoreMatrix);
	    std::swap(rowNames, colNames);
	    std::swap(gapExistenceCost, insExistenceCost);
	    std::swap(gapExtensionCost, insExtensionCost);
	  }
	  setLastSplitParams(params, opts,
			     scoreMatrix, rowNames.c_str(), colNames.c_str(),
			     gapExistenceCost, gapExtensionCost,
			     insExistenceCost, insExtensionCost,
			     scale, genomeSize, sequenceFormat);
	  opts.print();
	  params.print();
	  std::cout << "#\n";
	  state = 1;
	}
      }
      if (linePtr[0] == '#' && !startsWith(linePtr, "# batch")) {
	std::cout << linePtr << "\n";
      }
      if (state == 1) {  // we are reading alignments
	if (isBlankLine(linePtr)) {
	  addMaf(mafEnds, lineEnds);
	} else if (strchr(opts.no_split ? "asqpc" : "sqp", linePtr[0])) {
	  if (!opts.isTopSeqQuery && linePtr[0] == 's' && sLineCount++ % 2 &&
	      !isSameName(&inputText[qNameLineBeg], linePtr)) {
	    doOneBatch(inputText, lineEnds, mafEnds, splitter, opts, params,
		       isAlreadySplit);
	    eraseOldInput(inputText, lineEnds, mafEnds);
	    qNameLineBeg = lineEnds.back();
	  }
	  lineEnds.push_back(inputText.size());
	}
      }
      inputText.resize(lineEnds.back());
    }
  }
  addMaf(mafEnds, lineEnds);
  doOneBatch(inputText, lineEnds, mafEnds, splitter, opts, params,
	     isAlreadySplit);
}
