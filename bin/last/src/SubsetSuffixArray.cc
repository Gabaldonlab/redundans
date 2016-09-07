// Copyright 2008, 2009, 2010, 2013, 2014 Martin C. Frith

#include "SubsetSuffixArray.hh"
#include "io.hh"
#include <cassert>
#include <cstdio>  // remove
#include <sstream>

using namespace cbrc;

void SubsetSuffixArray::addPositions( const uchar* text,
				      indexT beg, indexT end, indexT step ){
  assert( step > 0 );
  const uchar* subsetMap = seed.firstMap();

  for( indexT i = beg; i < end; i += step ){
    if( subsetMap[ text[i] ] < CyclicSubsetSeed::DELIMITER ){
      suffixArray.v.push_back(i);
    }
    if( i + step < i ) break;  // avoid overflow
  }
}

void SubsetSuffixArray::clearPositions(){
  suffixArray.v.clear();
}

void SubsetSuffixArray::fromFiles( const std::string& baseName,
				   bool isMaskLowercase,
				   const uchar letterCode[] ){
  indexT textLength = 0;  // 0 never occurs in a valid file
  indexT unindexedPositions = 0;  // 0 never occurs in a valid file
  indexT bucketDepth = -1;
  seed.clear();

  std::string fileName = baseName + ".prj";
  std::ifstream f( fileName.c_str() );
  if( !f ) throw std::runtime_error( "can't open file: " + fileName );

  std::string line, word;
  while( getline( f, line ) ){
    std::istringstream iss(line);
    getline( iss, word, '=' );
    if( word == "totallength" ) iss >> textLength;
    if( word == "specialcharacters" ) iss >> unindexedPositions;
    if( word == "prefixlength" ) iss >> bucketDepth;
    if( word == "subsetseed" ){
      seed.appendPosition( iss, isMaskLowercase, letterCode );
    }
  }

  if( textLength == 0 || unindexedPositions == 0 || bucketDepth+1 == 0 ||
      !seed.span() || !f.eof() ){
    throw std::runtime_error( "can't read file: " + fileName );
  }

  indexT indexedPositions = textLength - unindexedPositions;
  suffixArray.m.open( baseName + ".suf", indexedPositions );
  makeBucketSteps( bucketDepth );
  buckets.m.open( baseName + ".bck", bucketSteps[0] );

  try{
    childTable.m.open( baseName + ".chi", indexedPositions );
  }catch( std::runtime_error ){
    try{
      kiddyTable.m.open( baseName + ".chi2", indexedPositions );
    }catch( std::runtime_error ){
      try{
	chibiTable.m.open( baseName + ".chi1", indexedPositions );
      }catch( std::runtime_error ){}
    }
  }
}

void SubsetSuffixArray::toFiles( const std::string& baseName,
				 bool isAppendPrj, indexT textLength ) const{
  assert( textLength > suffixArray.size() );

  std::string fileName = baseName + ".prj";
  std::ofstream f( fileName.c_str(),
		   isAppendPrj ? std::ios::app : std::ios::out );

  f << "totallength=" << textLength << '\n';
  f << "specialcharacters=" << textLength - suffixArray.size() << '\n';
  f << "prefixlength=" << maxBucketPrefix() << '\n';

  for( unsigned i = 0; i < seed.span(); ++i ){
    f << "subsetseed=";
    seed.writePosition( f, i );
    f << '\n';
  }

  f.close();
  if( !f ) throw std::runtime_error( "can't write file: " + fileName );

  memoryToBinaryFile( suffixArray.begin(), suffixArray.end(),
		      baseName + ".suf" );
  memoryToBinaryFile( buckets.begin(), buckets.end(), baseName + ".bck" );

  fileName = baseName + ".chi";
  std::remove( fileName.c_str() );
  memoryToBinaryFile( childTable.begin(), childTable.end(), fileName );

  fileName = baseName + ".chi2";
  std::remove( fileName.c_str() );
  memoryToBinaryFile( kiddyTable.begin(), kiddyTable.end(), fileName );

  fileName = baseName + ".chi1";
  std::remove( fileName.c_str() );
  memoryToBinaryFile( chibiTable.begin(), chibiTable.end(), fileName );
}

void SubsetSuffixArray::makeBuckets( const uchar* text, indexT bucketDepth ){
  if( bucketDepth+1 == 0 ) bucketDepth = defaultBucketDepth();

  makeBucketSteps( bucketDepth );

  buckets.v.clear();

  for( indexT i = 0; i < suffixArray.size(); ++i ){
    const uchar* textPtr = text + suffixArray[i];
    const uchar* subsetMap = seed.firstMap();
    indexT bucketIndex = 0;
    indexT depth = 0;

    while( depth < bucketDepth ){
      uchar subset = subsetMap[ *textPtr ];
      if( subset == CyclicSubsetSeed::DELIMITER ){
	bucketIndex += bucketSteps[depth] - 1;
	break;
      }
      ++textPtr;
      ++depth;
      indexT step = bucketSteps[depth];
      bucketIndex += subset * step;
      subsetMap = seed.nextMap( subsetMap );
    }

    buckets.v.resize( bucketIndex+1, i );
  }

  buckets.v.resize( bucketSteps[0], suffixArray.size() );
}

void SubsetSuffixArray::makeBucketSteps( indexT bucketDepth ){
  indexT step = 0;
  indexT depth = bucketDepth + 1;
  bucketSteps.resize( depth );

  while( depth > 0 ){
    --depth;
    step = step * seed.subsetCount(depth) + 1;
    bucketSteps[depth] = step;
  }
}

SubsetSuffixArray::indexT SubsetSuffixArray::defaultBucketDepth(){
  indexT maxBucketEntries = suffixArray.size() / 4;
  indexT bucketDepth = 0;
  indexT kmerEntries = 1;
  indexT bucketEntries = 1;

  while(true){
    indexT nextSubsetCount = seed.subsetCount(bucketDepth);
    if( kmerEntries > maxBucketEntries / nextSubsetCount ) return bucketDepth;
    kmerEntries *= nextSubsetCount;
    if( bucketEntries > maxBucketEntries - kmerEntries ) return bucketDepth;
    bucketEntries += kmerEntries;
    ++bucketDepth;
  }
}
