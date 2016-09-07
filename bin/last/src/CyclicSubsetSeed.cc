// Copyright 2009, 2010, 2013, 2014 Martin C. Frith

#include "CyclicSubsetSeed.hh"
#include "CyclicSubsetSeedData.hh"
#include "io.hh"
#include "stringify.hh"
#include <algorithm>  // sort
#include <sstream>
#include <stdexcept>
#include <cassert>
#include <cctype>  // toupper, tolower
#include <stddef.h>  // size_t
//#include <iostream>  // for debugging

#define ERR(x) throw std::runtime_error(x)

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

using namespace cbrc;

std::string CyclicSubsetSeed::stringFromName( const std::string& name ){
  for( std::size_t i = 0; i < COUNTOF(subsetSeeds); ++i )
    if( name == subsetSeeds[i].name )
      return subsetSeeds[i].text;

  return slurp( name );
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

bool CyclicSubsetSeed::nextPattern( std::istream& in,
				    std::vector< std::string >& seedAlphabet,
				    std::string& pattern ){
  while ( getline( in, pattern ) ){
    std::istringstream iss( pattern );
    std::string x, y;
    iss >> x;
    if( x.empty() || x[0] == '#' ) continue;
    iss >> y;
    if( y.empty() ) return true;
    if( x.size() > 1 ) ERR( "bad seed line: " + pattern );
    seedAlphabet.push_back( pattern );
  }
  return false;
}

static size_t findSeedLetter( const std::vector< std::string >& seedAlphabet,
			      char seedLetter ){
  // go backwards, so that newer definitions override older ones:
  for( size_t j = seedAlphabet.size(); j > 0; ){
    --j;
    const std::string& s = seedAlphabet[j];
    assert( !s.empty() );
    if( s[0] == seedLetter ) return j;
  }
  ERR( "unknown symbol in seed pattern: " + stringify(seedLetter) );
}

void CyclicSubsetSeed::init( const std::vector< std::string >& seedAlphabet,
			     const std::string& pattern,
			     bool isMaskLowercase,
			     const uchar letterCode[] ){
  clear();
  for( size_t i = 0; i < pattern.size(); ++i ){
    char seedLetter = pattern[i];
    if( std::isspace(seedLetter) ) continue;
    size_t j = findSeedLetter( seedAlphabet, seedLetter );
    std::istringstream iss( seedAlphabet[j] );
    iss >> seedLetter;
    appendPosition( iss, isMaskLowercase, letterCode );
  }
}

static void addLetter( uchar numbersToSubsets[], uchar letter, uchar subsetNum,
		       const uchar letterCode[] ){
  uchar number = letterCode[letter];
  if( number >= CyclicSubsetSeed::MAX_LETTERS )
    ERR( "bad symbol in subset-seed: " + stringify(letter) );
  if( numbersToSubsets[number] < CyclicSubsetSeed::DELIMITER )
    ERR( "repeated symbol in subset-seed: " + stringify(letter) );
  numbersToSubsets[number] = subsetNum;
}

void CyclicSubsetSeed::appendPosition( std::istream& inputLine,
				       bool isMaskLowercase,
				       const uchar letterCode[] ){
  std::string inputWord;
  std::vector<std::string> subsetList;
  std::vector<uchar> numbersToSubsets( MAX_LETTERS, DELIMITER );

  for( unsigned subsetNum = 0; inputLine >> inputWord; ++subsetNum ){
    assert( subsetNum < DELIMITER );
    std::string subset;

    for( size_t i = 0; i < inputWord.size(); ++i ){
      uchar upper = std::toupper( inputWord[i] );
      uchar lower = std::tolower( inputWord[i] );
      addLetter( &numbersToSubsets[0], upper, subsetNum, letterCode );
      subset += upper;
      if( !isMaskLowercase && lower != upper ){
	addLetter( &numbersToSubsets[0], lower, subsetNum, letterCode );
      }
    }

    sort( subset.begin(), subset.end() );  // canonicalize
    subsetList.push_back( subset );
  }

  subsetLists.push_back( subsetList );
  subsetMaps.insert( subsetMaps.end(),
		     numbersToSubsets.begin(), numbersToSubsets.end() );
}

void CyclicSubsetSeed::writePosition( std::ostream& out,
				      unsigned position ) const{
  assert( position < subsetLists.size() );
  for( size_t i = 0; i < subsetLists[position].size(); ++i ){
    if( i > 0 ) out << ' ';
    out << subsetLists[position][i];
  }
}
