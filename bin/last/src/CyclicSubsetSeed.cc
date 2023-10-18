// Copyright 2009, 2010, 2013, 2014 Martin C. Frith

#include "CyclicSubsetSeed.hh"
#include "CyclicSubsetSeedData.hh"
#include "zio.hh"
#include "stringify.hh"
#include <algorithm>  // sort
#include <sstream>
#include <cassert>
#include <cctype>  // toupper, tolower
#include <string.h>
//#include <iostream>  // for debugging

#define ERR(x) throw std::runtime_error(x)

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

using namespace cbrc;

static const char *canonicalName(const char *name) {
  for (size_t i = 0; i < COUNTOF(subsetSeedNicknames); ++i)
    if (strcmp(name, subsetSeedNicknames[i].nickname) == 0)
      return subsetSeedNicknames[i].realname;
  return name;
}

std::string CyclicSubsetSeed::stringFromName( const std::string& name ){
  const char *n = canonicalName(name.c_str());

  for (size_t i = 0; i < COUNTOF(subsetSeeds); ++i)
    if (strcmp(n, subsetSeeds[i].name) == 0)
      return subsetSeeds[i].text;

  return slurp(n);
}

std::string
CyclicSubsetSeed::stringFromPatterns( const std::string& patterns,
				      const std::string& sequenceLetters ){
  std::string spacedLetters;
  for( size_t i = 0; i < sequenceLetters.size(); ++i ){
    if( i > 0 ) spacedLetters += ' ';
    spacedLetters += sequenceLetters[i];
  }

  std::string p = patterns;
  for( size_t i = 0; i < p.size(); ++i )
    switch( p[i] ){
    case ',':
      p[i] = '\n'; break;
    case '#':
      p[i] = '1'; break;
    case '_':
    case '-':
      p[i] = '0'; break;
    case 't':
    case '@':
      p[i] = 'T'; break;
    }

  return "\
1  " + spacedLetters + "\n\
0  " + sequenceLetters + "\n\
T  AG CT\n\
" + p;
}

std::string
CyclicSubsetSeed::stringFromDnaPatterns(std::string patterns) {
  for (size_t i = 0; i < patterns.size(); ++i) {
    if (patterns[i] == ',') patterns[i] = '\n';
  }

  return "\
A  A\n\
C  C\n\
G  G\n\
T  T\n\
R  A G\n\
r  AG\n\
Y  C T\n\
y  CT\n\
N  A C G T\n\
n  ACGT\n\
@  AG CT\n\
" + patterns;
}

static int lineType(const std::string &line) {
  std::istringstream s(line);
  std::string x, y;
  s >> x >> y;
  return (x.empty() || x[0] == '#')    ? 0  // blank or comment line
    :    (x.size() == 1 && !y.empty()) ? 1  // seed-alphabet line
    : 2;                                    // seed-pattern line
}

static const char *letterGroups(const std::vector<std::string> &seedAlphabet,
				char seedLetter) {
  // go backwards, so that newer definitions override older ones:
  for (size_t i = seedAlphabet.size(); i-- > 0; ) {
    const char *s = seedAlphabet[i].c_str();
    if (s[0] == seedLetter) return s + 2;
  }
  ERR("unknown symbol in seed pattern: " + stringify(seedLetter));
}

void CyclicSubsetSeed::addPatterns(std::vector<CyclicSubsetSeed> &patterns,
				   const std::string &text,
				   bool isMaskLowercase,
				   const uchar letterCode[],
				   const std::string &mainSequenceAlphabet) {
  std::vector<std::string> seedAlphabet;
  std::string line, word;
  std::istringstream textStream(text);

  while (getline(textStream, line)) {
    int type = lineType(line);
    if (type == 1) {
      seedAlphabet.push_back(line);
    } else if (type == 2) {
      std::istringstream lineStream(line);
      while (lineStream >> word) {
	CyclicSubsetSeed pat;
	for (size_t i = 0; i < word.size(); ++i) {
	  std::istringstream s(letterGroups(seedAlphabet, word[i]));
	  pat.appendPosition(s, isMaskLowercase, letterCode,
			     mainSequenceAlphabet);
	}
	patterns.push_back(pat);
      }
    }
  }
}

static void addLetter( uchar toSubsetNum[], uchar letter, unsigned subsetNum,
		       const uchar letterCode[] ){
  uchar number = letterCode[letter];
  if( number >= CyclicSubsetSeed::MAX_LETTERS )
    ERR( "bad symbol in subset-seed: " + stringify(letter) );
  if( toSubsetNum[number] < CyclicSubsetSeed::DELIMITER )
    ERR( "repeated symbol in subset-seed: " + stringify(letter) );
  toSubsetNum[number] = subsetNum;
}

static void addLowercase(bool isMaskLowercase, uchar *toSubsetNum, uchar upper,
			 unsigned subsetNum, const uchar *letterCode) {
  uchar lower = std::tolower(upper);
  if (!isMaskLowercase && lower != upper) {
    addLetter(toSubsetNum, lower, subsetNum, letterCode);
  }
}

void CyclicSubsetSeed::appendPosition(std::istream& inputLine,
				      bool isMaskLowercase,
				      const uchar letterCode[],
				      const std::string &mainSequenceAlphabet){
  std::string inputWord;
  std::vector<std::string> subsetList;
  std::vector<uchar> toSubsetNum(MAX_LETTERS, DELIMITER);
  unsigned subsetNum = 0;
  unsigned letterNum = 0;

  while (inputLine >> inputWord) {
    assert( subsetNum < DELIMITER );
    std::string subset;

    for( size_t i = 0; i < inputWord.size(); ++i ){
      uchar c = inputWord[i];
      uchar u = std::toupper(c);
      addLetter(&toSubsetNum[0], u, subsetNum, letterCode);
      addLowercase(isMaskLowercase, &toSubsetNum[0], u, subsetNum, letterCode);
      ++letterNum;
      subset += u;
    }

    sort( subset.begin(), subset.end() );  // canonicalize
    subsetList.push_back( subset );
    ++subsetNum;
  }

  const bool isExactSeed = (letterNum == subsetNum);

  bool isNewSubset = false;
  for (size_t i = 0; i < mainSequenceAlphabet.size(); ++i) {
    uchar c = mainSequenceAlphabet[i];
    uchar u = std::toupper(c);
    uchar number = letterCode[u];
    assert(number < MAX_LETTERS);
    if (toSubsetNum[number] == DELIMITER) {
      toSubsetNum[number] = subsetNum;
      addLowercase(isMaskLowercase, &toSubsetNum[0], u, subsetNum, letterCode);
      subsetNum += isExactSeed;
      isNewSubset = !isExactSeed;
    }
  }
  subsetNum += isNewSubset;

  subsetLists.push_back( subsetList );
  subsetMaps.insert(subsetMaps.end(), toSubsetNum.begin(), toSubsetNum.end());
  originalSubsetMaps.insert(originalSubsetMaps.end(),
			    toSubsetNum.begin(), toSubsetNum.end());
  numOfSubsetsPerPosition.push_back(subsetNum);
}

void CyclicSubsetSeed::writePosition( std::ostream& out,
				      size_t position ) const{
  assert( position < subsetLists.size() );
  for( size_t i = 0; i < subsetLists[position].size(); ++i ){
    if( i > 0 ) out << ' ';
    out << subsetLists[position][i];
  }
}
