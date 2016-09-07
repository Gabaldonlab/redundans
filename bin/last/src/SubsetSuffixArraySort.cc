// Copyright 2008, 2009, 2010, 2011, 2013 Martin C. Frith

// Parts of this code are adapted from "Engineering Radix Sort" by PM
// McIlroy, K Bostic, MD McIlroy.

#include "SubsetSuffixArray.hh"
#include <algorithm>  // iter_swap, min
//#include <iostream>  // for debugging

using namespace cbrc;

namespace{
  typedef SubsetSuffixArray::indexT indexT;
  struct Stack{ indexT* beg; indexT* end; indexT depth; };
  Stack stack[1048576];  // big enough???
  Stack* sp = stack;
}

#define PUSH(b, e, d) if( e-b > 1 ) sp->beg = b, sp->end = e, (sp++)->depth = d
#define  POP(b, e, d) b = (--sp)->beg, e = sp->end, d = sp->depth

static void insertionSort( const uchar* text, const CyclicSubsetSeed& seed,
			   indexT* beg, indexT* end, const uchar* subsetMap ){
  for( indexT* i = beg+1; i < end; ++i ){
    const uchar* newText = text + *i;
    for( indexT* j = i; j > beg; --j ){
      indexT* k = j - 1;
      const uchar* oldText = text + *k;
      if( !seed.isLess( newText, oldText, subsetMap ) ) break;
      std::iter_swap( j, k );
    }
  }
}

void SubsetSuffixArray::sort2( const uchar* text, indexT* beg,
			       const uchar* subsetMap ){
  indexT* mid = beg + 1;

  const uchar* s = text + *beg;
  const uchar* t = text + *mid;
  while( true ){
    uchar x = subsetMap[ *s ];
    uchar y = subsetMap[ *t ];
    if( x != y ){
      if( x > y ) std::iter_swap( beg, mid );
      break;
    }
    if( x == CyclicSubsetSeed::DELIMITER ) return;
    ++s;
    ++t;
    subsetMap = seed.nextMap(subsetMap);
  }

  if( isChildDirectionForward( beg ) ){
    setChildForward( beg + 0, mid );
  }else{
    setChildReverse( beg + 2, mid );
  }
}

// Specialized sort for 1 symbol + 1 delimiter.
// E.g. wildcard positions in spaced seeds.
void SubsetSuffixArray::radixSort1( const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* begN = end;  // beginning of delimiters

  while( end0 < begN ){
    const indexT x = *end0;
    switch( subsetMap[ text[x] ] ){
    case 0:
      end0++;
      break;
    default:  // the delimiter subset
      *end0 = *--begN;
      *begN = x;
      break;
    }
  }

  PUSH( beg, end0, depth );   // the '0's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
  }
}

// Specialized sort for 2 symbols + 1 delimiter.
// E.g. transition-constrained positions in subset seeds.
void SubsetSuffixArray::radixSort2( const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* begN = end;  // beginning of delimiters

  while( end1 < begN ){
    const indexT x = *end1;
    switch( subsetMap[ text[x] ] ){
      case 0:
        *end1++ = *end0;
        *end0++ = x;
        break;
      case 1:
        end1++;
        break;
      default:  // the delimiter subset
        *end1 = *--begN;
        *begN = x;
        break;
    }
  }

  PUSH( beg, end0, depth );   // the '0's
  PUSH( end0, end1, depth );  // the '1's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

// Specialized sort for 3 symbols + 1 delimiter.
// E.g. subset seeds for bisulfite-converted DNA.
void SubsetSuffixArray::radixSort3( const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* beg2 = end;  // beginning of '2's
  indexT* begN = end;  // beginning of delimiters

  while( end1 < beg2 ){
    const indexT x = *end1;
    switch( subsetMap[ text[x] ] ){
      case 0:
        *end1++ = *end0;
        *end0++ = x;
        break;
      case 1:
        end1++;
        break;
      case 2:
        *end1 = *--beg2;
        *beg2 = x;
        break;
      default:  // the delimiter subset
        *end1 = *--beg2;
        *beg2 = *--begN;
        *begN = x;
        break;
    }
  }

  PUSH( beg, end0, depth );   // the '0's
  PUSH( end0, end1, depth );  // the '1's
  PUSH( beg2, begN, depth );  // the '2's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
    if( begN == end ) return;
    setChildForward( beg2, begN );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( beg2 == beg ) return;
    setChildReverse( begN, beg2 );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

// Specialized sort for 4 symbols + 1 delimiter.  E.g. DNA.
void SubsetSuffixArray::radixSort4( const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth ){
  indexT* end0 = beg;  // end of '0's
  indexT* end1 = beg;  // end of '1's
  indexT* end2 = beg;  // end of '2's
  indexT* beg3 = end;  // beginning of '3's
  indexT* begN = end;  // beginning of delimiters

  while( end2 < beg3 ){
    const indexT x = *end2;
    switch( subsetMap[ text[x] ] ){
    case 0:
      *end2++ = *end1;
      *end1++ = *end0;
      *end0++ = x;
      break;
    case 1:
      *end2++ = *end1;
      *end1++ = x;
      break;
    case 2:
      end2++;
      break;
    case 3:
      *end2 = *--beg3;
      *beg3 = x;
      break;
    default:  // the delimiter subset
      *end2 = *--beg3;
      *beg3 = *--begN;
      *begN = x;
      break;
    }
  }

  PUSH( beg, end0, depth );   // the '0's
  PUSH( end0, end1, depth );  // the '1's
  PUSH( end1, end2, depth );  // the '2's
  PUSH( beg3, begN, depth );  // the '3's

  if( isChildDirectionForward( beg ) ){
    if( end0 == end ) return;
    setChildForward( beg, end0 );
    if( end1 == end ) return;
    setChildForward( end0, end1 );
    if( end2 == end ) return;
    setChildForward( end1, end2 );
    if( begN == end ) return;
    setChildForward( beg3, begN );
  }else{
    if( begN == beg ) return;
    setChildReverse( end, begN );
    if( beg3 == beg ) return;
    setChildReverse( begN, beg3 );
    if( end1 == beg ) return;
    setChildReverse( end2, end1 );
    if( end0 == beg ) return;
    setChildReverse( end1, end0 );
  }
}

void SubsetSuffixArray::radixSortN( const uchar* text, const uchar* subsetMap,
				    indexT* beg, indexT* end, indexT depth,
				    unsigned subsetCount ){
  static indexT bucketSize[256];  // initialized to zero at startup
  /*  */ indexT* bucketEnd[256];  // "static" makes little difference to speed

  // get bucket sizes (i.e. letter counts):
  // The intermediate oracle array makes it faster (see "Engineering
  // Radix Sort for Strings" by J Karkkainen & T Rantala)
  for( indexT* i = beg; i < end; /* noop */ ){
    uchar oracle[256];
    uchar* oracleEnd =
      oracle + std::min( sizeof(oracle), std::size_t(end - i) );
    for( uchar* j = oracle; j < oracleEnd; ++j )
      *j = subsetMap[ text[ *i++ ] ];
    for( uchar* j = oracle; j < oracleEnd; ++j )
      ++bucketSize[ *j ];
  }

  // get bucket ends, and put buckets on the stack to sort within them later:
  // (could push biggest bucket first, to ensure logarithmic stack growth)
  indexT* pos = beg;
  for( unsigned i = 0; i < subsetCount; ++i ){
    indexT* nextPos = pos + bucketSize[i];
    PUSH( pos, nextPos, depth );
    pos = nextPos;
    bucketEnd[i] = pos;
  }
  // don't sort within the delimiter bucket:
  bucketEnd[ CyclicSubsetSeed::DELIMITER ] = end;

  if( isChildDirectionForward( beg ) ){
    pos = beg;
    for( unsigned i = 0; i < subsetCount; ++i ){
      indexT* nextPos = bucketEnd[i];
      if( nextPos == end ) break;
      setChildForward( pos, nextPos );
      pos = nextPos;
    }
  }else{
    pos = end;
    for( unsigned i = subsetCount; i > 0; --i ){
      indexT* nextPos = bucketEnd[i - 1];
      if( nextPos == beg ) break;
      setChildReverse( pos, nextPos );
      pos = nextPos;
    }
  }

  // permute items into the correct buckets:
  for( indexT* i = beg; i < end; /* noop */ ) {
    unsigned subset;  // unsigned is faster than uchar!
    indexT holdOut = *i;
    while( --bucketEnd[ subset = subsetMap[ text[holdOut] ] ] > i ){
      std::swap( *bucketEnd[subset], holdOut );
    }
    *i = holdOut;
    i += bucketSize[subset];
    bucketSize[subset] = 0;  // reset it so we can reuse it
  }
}

void SubsetSuffixArray::sortIndex( const uchar* text,
				   indexT maxUnsortedInterval,
				   int childTableType ){
  const indexT minLength = 1;

  if( childTableType == 1 ) chibiTable.v.assign( suffixArray.v.size(), -1 );
  if( childTableType == 2 ) kiddyTable.v.assign( suffixArray.v.size(), -1 );
  if( childTableType == 3 ) childTable.v.assign( suffixArray.v.size(), 0 );

  PUSH( &suffixArray.v.front(), &suffixArray.v.back() + 1, 0 );
  setChildReverse( &suffixArray.v.back() + 1, &suffixArray.v.front() );

  while( sp > stack ){
    indexT* beg;
    indexT* end;
    indexT depth;
    POP( beg, end, depth );

    if( end - beg <= maxUnsortedInterval && depth >= minLength ) continue;

    const uchar* textBase = text + depth;
    const uchar* subsetMap = seed.subsetMap(depth);

    if( childTableType == 0 ){
      if( end - beg < 10 ){  // ???
	insertionSort( textBase, seed, beg, end, subsetMap );
	continue;
      }
    }else{
      if( end - beg == 2 ){
	sort2( textBase, beg, subsetMap );
	continue;
      }
    }

    unsigned subsetCount = seed.subsetCount(depth);

    ++depth;

    switch( subsetCount ){
      case 1:  radixSort1( textBase, subsetMap, beg, end, depth );  break;
      case 2:  radixSort2( textBase, subsetMap, beg, end, depth );  break;
      case 3:  radixSort3( textBase, subsetMap, beg, end, depth );  break;
      case 4:  radixSort4( textBase, subsetMap, beg, end, depth );  break;
      default: radixSortN( textBase, subsetMap, beg, end, depth, subsetCount );
    }
  }
}
