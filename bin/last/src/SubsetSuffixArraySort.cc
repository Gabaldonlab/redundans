// Copyright 2008, 2009, 2010, 2011, 2013 Martin C. Frith

// Parts of this code are adapted from "Engineering Radix Sort" by PM
// McIlroy, K Bostic, MD McIlroy.

#include "SubsetSuffixArray.hh"
#include <algorithm>  // swap, min
//#include <iostream>  // for debugging

#ifdef HAS_CXX_THREADS
#include <thread>
#endif

using namespace cbrc;

static void posCpy(PosPart *dest, PosPart *src) {
  for (int i = 0; i < posParts; ++i) dest[i] = src[i];
}

static void posSwap(PosPart *x, PosPart *y) {
  for (int i = 0; i < posParts; ++i) std::swap(x[i], y[i]);
}

namespace{
  typedef SubsetSuffixArray::indexT indexT;
  typedef SubsetSuffixArray::Range Range;
}

static void pushRange(std::vector<Range> &v,
		      PosPart *beg, PosPart *end, indexT depth) {
  if (end - beg > posParts) {
    Range r = {beg, end, depth};
    v.push_back(r);
  }
}

static void insertionSort(const uchar *text, const CyclicSubsetSeed &seed,
			  PosPart *beg, PosPart *end, const uchar *subsetMap) {
  for (PosPart *i = beg + posParts; i < end; i += posParts) {
    const uchar *newText = text + posGet(i);
    PosPart *j = i;
    do {
      PosPart *k = j;
      j -= posParts;
      const uchar *oldText = text + posGet(j);
      if( !seed.isLess( newText, oldText, subsetMap ) ) break;
      posSwap(j, k);
    } while (j > beg);
  }
}

void SubsetSuffixArray::sort2(const uchar* text, const CyclicSubsetSeed &seed,
			      PosPart *beg, const uchar* subsetMap) {
  PosPart *mid = beg + posParts;

  const uchar *s = text + posGet(beg);
  const uchar *t = text + posGet(mid);
  while( true ){
    uchar x = subsetMap[ *s ];
    uchar y = subsetMap[ *t ];
    if( x != y ){
      if (x > y) posSwap(beg, mid);
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
    setChildReverse( beg + 2 * posParts, mid );
  }
}

// Specialized sort for 1 symbol + 1 delimiter.
// E.g. wildcard positions in spaced seeds.
void SubsetSuffixArray::radixSort1( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    PosPart *beg, PosPart *end, indexT depth ){
  PosPart *end0 = beg;  // end of '0's
  PosPart *begN = end;  // beginning of delimiters

  while( end0 < begN ){
    const size_t x = posGet(end0);
    switch( subsetMap[ text[x] ] ){
    case 0:
      end0 += posParts;
      break;
    default:  // the delimiter subset
      begN -= posParts;  posCpy(end0, begN);
      posSet(begN, x);
      break;
    }
  }

  pushRange( rangeStack, beg, end0, depth );   // the '0's

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
void SubsetSuffixArray::radixSort2( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    PosPart *beg, PosPart *end, indexT depth ){
  PosPart *end0 = beg;  // end of '0's
  PosPart *end1 = beg;  // end of '1's
  PosPart *begN = end;  // beginning of delimiters

  while( end1 < begN ){
    const size_t x = posGet(end1);
    switch( subsetMap[ text[x] ] ){
      case 0:
        posCpy(end1, end0);  end1 += posParts;
        posSet(end0, x);     end0 += posParts;
        break;
      case 1:
        end1 += posParts;
        break;
      default:  // the delimiter subset
        begN -= posParts;  posCpy(end1, begN);
	posSet(begN, x);
        break;
    }
  }

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's

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
void SubsetSuffixArray::radixSort3( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    PosPart *beg, PosPart *end, indexT depth ){
  PosPart *end0 = beg;  // end of '0's
  PosPart *end1 = beg;  // end of '1's
  PosPart *beg2 = end;  // beginning of '2's
  PosPart *begN = end;  // beginning of delimiters

  while( end1 < beg2 ){
    const size_t x = posGet(end1);
    switch( subsetMap[ text[x] ] ){
      case 0:
        posCpy(end1, end0);  end1 += posParts;
        posSet(end0, x);     end0 += posParts;
        break;
      case 1:
	end1 += posParts;
        break;
      case 2:
        beg2 -= posParts;  posCpy(end1, beg2);
        posSet(beg2, x);
        break;
      default:  // the delimiter subset
        beg2 -= posParts;  posCpy(end1, beg2);
        begN -= posParts;  posCpy(beg2, begN);
        posSet(begN, x);
        break;
    }
  }

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's
  pushRange( rangeStack, beg2, begN, depth );  // the '2's

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
void SubsetSuffixArray::radixSort4( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    PosPart *beg, PosPart *end, indexT depth ){
  PosPart *end0 = beg;  // end of '0's
  PosPart *end1 = beg;  // end of '1's
  PosPart *end2 = beg;  // end of '2's
  PosPart *beg3 = end;  // beginning of '3's
  PosPart *begN = end;  // beginning of delimiters

  while( end2 < beg3 ){
    const size_t x = posGet(end2);
    switch( subsetMap[ text[x] ] ){
    case 0:
      posCpy(end2, end1);  end2 += posParts;
      posCpy(end1, end0);  end1 += posParts;
      posSet(end0, x);     end0 += posParts;
      break;
    case 1:
      posCpy(end2, end1);  end2 += posParts;
      posSet(end1, x);     end1 += posParts;
      break;
    case 2:
      end2 += posParts;
      break;
    case 3:
      beg3 -= posParts;  posCpy(end2, beg3);
      posSet(beg3, x);
      break;
    default:  // the delimiter subset
      beg3 -= posParts;  posCpy(end2, beg3);
      begN -= posParts;  posCpy(beg3, begN);
      posSet(begN, x);
      break;
    }
  }

  pushRange( rangeStack, beg, end0, depth );   // the '0's
  pushRange( rangeStack, end0, end1, depth );  // the '1's
  pushRange( rangeStack, end1, end2, depth );  // the '2's
  pushRange( rangeStack, beg3, begN, depth );  // the '3's

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

const unsigned numOfBuckets = 256;

void SubsetSuffixArray::radixSortN( std::vector<Range>& rangeStack,
				    const uchar* text, const uchar* subsetMap,
				    PosPart *beg, PosPart *end, indexT depth,
				    unsigned subsetCount, indexT* bucketSize ){
  PosPart *bucketEnd[numOfBuckets];

  // get bucket sizes (i.e. letter counts):
  // The intermediate oracle array makes it faster (see "Engineering
  // Radix Sort for Strings" by J Karkkainen & T Rantala)
  for (PosPart *i = beg; i < end; /* noop */) {
    uchar oracle[256];
    PosPart *iEnd = i + std::min(sizeof(oracle) * posParts, size_t(end - i));
    uchar *j = oracle;
    while (i < iEnd) {
      *j++ = subsetMap[text[posGet(i)]];
      i += posParts;
    }
    for (uchar *k = oracle; k < j; ++k) {
      bucketSize[*k] += posParts;
    }
  }

  // get bucket ends, and put buckets on the stack to sort within them later:
  // (could push biggest bucket first, to ensure logarithmic stack growth)
  PosPart *pos = beg;
  for( unsigned i = 0; i < subsetCount; ++i ){
    PosPart *nextPos = pos + bucketSize[i];
    pushRange( rangeStack, pos, nextPos, depth );
    pos = nextPos;
    bucketEnd[i] = pos;
  }
  // don't sort within the delimiter bucket:
  bucketEnd[ CyclicSubsetSeed::DELIMITER ] = end;

  if( isChildDirectionForward( beg ) ){
    pos = beg;
    for( unsigned i = 0; i < subsetCount; ++i ){
      PosPart *nextPos = bucketEnd[i];
      if( nextPos == end ) break;
      setChildForward( pos, nextPos );
      pos = nextPos;
    }
  }else{
    pos = end;
    for( unsigned i = subsetCount; i > 0; --i ){
      PosPart *nextPos = bucketEnd[i - 1];
      if( nextPos == beg ) break;
      setChildReverse( pos, nextPos );
      pos = nextPos;
    }
  }

  // permute items into the correct buckets:
  for (PosPart *i = beg; i < end; /* noop */) {
    unsigned subset;  // unsigned is faster than uchar!
    while (1) {
      subset = subsetMap[text[posGet(i)]];
      bucketEnd[subset] -= posParts;
      if (bucketEnd[subset] <= i) break;
      posSwap(i, bucketEnd[subset]);
    }
    i += bucketSize[subset];
    bucketSize[subset] = 0;  // reset it so we can reuse it
  }
}

static size_t rangeSize(const Range &r) {
  return r.end - r.beg;
}

static size_t nextRangeSize(const std::vector<Range> &ranges) {
  return rangeSize(ranges.back());
}

static size_t rangeSizeSum(const std::vector<Range> &ranges) {
  size_t s = 0;
  for (size_t i = 0; i < ranges.size(); ++i) {
    s += rangeSize(ranges[i]);
  }
  return s;
}

static size_t numOfThreadsForOneRange(size_t numOfThreads,
				      size_t sizeOfThisRange,
				      size_t sizeOfAllRanges) {
  // We want:
  // min(t | sizeOfThisRange / t < sizeOfOtherRanges / (numOfThreads - (t+1)))
  // Or equivalently:
  // max(t | sizeOfThisRange / (t-1) >= sizeOfOtherRanges / (numOfThreads - t))
  double x = numOfThreads - 1;  // use double to avoid overflow
  return (x * sizeOfThisRange + sizeOfAllRanges) / sizeOfAllRanges;
}

void SubsetSuffixArray::sortRanges(std::vector<Range> *stacks,
				   indexT *bucketSizes,
				   const uchar *text,
				   unsigned wordLength,
				   const CyclicSubsetSeed &seed,
				   size_t maxUnsortedInterval,
				   size_t numOfThreads) {
  std::vector<Range> &myStack = stacks[0];

  while( !myStack.empty() ){
#ifdef HAS_CXX_THREADS
    size_t numOfChunks = std::min(numOfThreads, myStack.size());
    if (numOfChunks > 1) {
      size_t totalSize = rangeSizeSum(myStack);
      size_t numOfNewThreads = numOfChunks - 1;
      std::vector<std::thread> threads(numOfNewThreads);

      for (size_t i = 0; i < numOfNewThreads; ++i) {
	size_t thisSize = nextRangeSize(myStack);
	size_t t = numOfThreadsForOneRange(numOfThreads, thisSize, totalSize);
	size_t maxThreads = numOfThreads - (numOfNewThreads - i);
	size_t thisThreads = std::min(t, maxThreads);
	numOfThreads -= thisThreads;

	do {
	  totalSize -= nextRangeSize(myStack);
	  stacks[numOfThreads].push_back(myStack.back());
	  myStack.pop_back();
	  thisSize += nextRangeSize(myStack);
	} while (myStack.size() > numOfThreads &&
		 thisSize <= totalSize / numOfThreads);
	// We want:
	// max(r | sizeSum(r) <= (totalSize - sizeSum(r-1)) / newNumOfThreads)

	threads[i] = std::thread(&SubsetSuffixArray::sortRanges, this,
				 stacks + numOfThreads,
				 bucketSizes + numOfThreads * numOfBuckets,
				 text, wordLength, seed,
				 maxUnsortedInterval, thisThreads);
      }
      sortRanges(stacks, bucketSizes, text, wordLength, seed,
		 maxUnsortedInterval, numOfThreads);
      for (size_t i = 0; i < numOfNewThreads; ++i) {
	threads[i].join();
      }
      return;
    }
#endif

    PosPart *beg = myStack.back().beg;
    PosPart *end = myStack.back().end;
    indexT depth = myStack.back().depth;
    myStack.pop_back();

    size_t interval = posCount(beg, end);
    const indexT minLength = 1;
    if( interval <= maxUnsortedInterval && depth >= minLength ) continue;

    const uchar* textBase = text + depth;
    const uchar* subsetMap = seed.subsetMap(depth);

    if( childTable.v.empty() && kiddyTable.v.empty() && chibiTable.v.empty() ){
      if( interval < 10 ){  // ???
	insertionSort( textBase, seed, beg, end, subsetMap );
	continue;
      }
    }else{
      if( interval == 2 ){
	sort2( textBase, seed, beg, subsetMap );
	continue;
      }
    }

    unsigned subsetCount;
    if (depth < wordLength) {
      subsetCount = seed.restrictedSubsetCount(depth);
      // xxx inefficient
    } else {
      subsetCount = seed.unrestrictedSubsetCount(depth);
    }

    ++depth;

    switch( subsetCount ){
    case 1:  radixSort1(myStack, textBase, subsetMap, beg, end, depth); break;
    case 2:  radixSort2(myStack, textBase, subsetMap, beg, end, depth); break;
    case 3:  radixSort3(myStack, textBase, subsetMap, beg, end, depth); break;
    case 4:  radixSort4(myStack, textBase, subsetMap, beg, end, depth); break;
    default: radixSortN(myStack, textBase, subsetMap, beg, end, depth,
			subsetCount, bucketSizes);
    }
  }
}

void SubsetSuffixArray::sortIndex( const uchar* text,
				   unsigned wordLength,
				   const size_t* cumulativeCounts,
				   size_t maxUnsortedInterval,
				   int childTableType,
				   size_t numOfThreads ){
  if (childTableType == 1) chibiTable.v.assign(size(), -1);
  if (childTableType == 2) kiddyTable.v.assign(size(), -1);
  if (childTableType == 3) childTable.v.assign(size(), 0);

  std::vector< std::vector<Range> > stacks(numOfThreads);
  std::vector<indexT> bucketSizes(numOfThreads * numOfBuckets);
  PosPart *a = &suffixArray.v[0];

  PosPart *beg = a;
  for (size_t i = 0; i < seeds.size(); ++i) {
    PosPart *end = a + cumulativeCounts[i] * posParts;
    pushRange(stacks[0], beg, end, 0);
    setChildReverse(end, beg);
    sortRanges(&stacks[0], &bucketSizes[0], text, wordLength, seeds[i],
	       maxUnsortedInterval, numOfThreads);
    beg = end;
  }
}
