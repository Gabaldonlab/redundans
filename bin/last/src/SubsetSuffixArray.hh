// Copyright 2008, 2009, 2010, 2013, 2014 Martin C. Frith

// This class holds a suffix array.  The suffix array is just a list
// of numbers indicating positions in a text, sorted according to the
// alphabetical order of the text suffixes starting at these
// positions.  A query sequence can then be matched incrementally to
// the suffix array using binary search.

// A "subset suffix array" means that, when comparing two suffixes, we
// consider subsets of letters to be equivalent.  For example, we
// might consider purines to be equivalent to each other, and
// pyrimidines to be equivalent to each other.  The subsets may vary
// from position to position as we compare two suffixes.

// There is always a special subset, called DELIMITER, which doesn't
// match anything.

// For faster matching, we use "buckets", which store the start and
// end in the suffix array of all size-k prefixes of the suffixes.
// They store this information for all values of k from 1 to, say, 12.

// This class can store multiple concatenated suffix arrays: each
// array holds suffixes starting with each pattern (of length
// "wordLength") in a DnaWordsFinder.  The endpoints of the
// concatenated arrays are specified by "cumulativeCounts".  Each
// array has its own letter-subsets.  "seedNum" specifies one of the
// arrays.

#ifndef SUBSET_SUFFIX_ARRAY_HH
#define SUBSET_SUFFIX_ARRAY_HH

#include "CyclicSubsetSeed.hh"
#include "dna_words_finder.hh"
#include "VectorOrMmap.hh"

#include <climits>

namespace cbrc{

#if LAST_POS_BYTES == 8
  typedef size_t PosPart;
  const int posParts = 1;
#elif LAST_POS_BYTES == 5
  typedef unsigned char PosPart;
  const int posParts = 5;
#else
  typedef unsigned PosPart;
  const int posParts = 1;
#endif

typedef PosPart OffPart;
const int offParts = posParts;

inline size_t posGet(const PosPart *p) {
  size_t x = 0;
  for (int i = 0; i < posParts; ++i) {
    size_t y = p[i];  // must convert to size_t before shifting!
    x += y << (i * sizeof(PosPart) * CHAR_BIT);
  }
  return x;
}

inline void posSet(PosPart *p, size_t value) {
  for (int i = 0; i < posParts; ++i) {
    p[i] = value >> (i * sizeof(PosPart) * CHAR_BIT);
  }
}

inline size_t posCount(const PosPart *beg, const PosPart *end) {
  size_t d = end - beg;
  return d / posParts;  // faster if the dividend is unsigned?
}

class SubsetSuffixArray{
public:

#if LAST_POS_BYTES > 4
  typedef size_t indexT;
#else
  typedef unsigned indexT;
#endif

  struct Range {PosPart *beg; PosPart *end; indexT depth;};

  std::vector<CyclicSubsetSeed> &getSeeds() { return seeds; }
  const std::vector<CyclicSubsetSeed> &getSeeds() const { return seeds; }

  size_t size() const { return suffixArray.size() / posParts; }

  PosPart *resizedPositions(size_t numOfPositions) {
    suffixArray.v.resize(numOfPositions * posParts);
    return &suffixArray.v[0];
  }

  // Add positions in the range [beg,end) that are "minimizers" for
  // the given window and seed pattern.  (Only minimizers at each
  // step-th position are added: this step parameter may be useless.)
  // Positions starting with delimiters aren't added.
  // The positions aren't sorted.
  void addMinimizerPositions(const uchar *text, size_t beg, size_t end,
			     size_t step, size_t minimizerWindow);

  // Store positions in [seqBeg, seqEnd) where certain "words" start.
  // The cumulative word counts must be provided.  (cumulativeCounts
  // is internally modified and restored to its original values).
  void setWordPositions(const DnaWordsFinder &finder, size_t *cumulativeCounts,
			const uchar *seqBeg, const uchar *seqEnd);

  // Sort the suffix array (but don't make the buckets).
  void sortIndex(const uchar *text,
		 unsigned wordLength, const size_t *cumulativeCounts,
		 size_t maxUnsortedInterval, int childTableType,
		 size_t numOfThreads);

  // Make the buckets.  If bucketDepth+1 == 0, then the bucket depth
  // is: the maximum possible such that (memory use of buckets) <=
  // (memory use of stored positions) / minPositionsPerBucket.
  void makeBuckets(const uchar *text,
		   unsigned wordLength, const size_t *cumulativeCounts,
		   size_t minPositionsPerBucket, unsigned bucketDepth,
		   size_t numOfThreads);

  void fromFiles(const std::string &baseName,
		 bool isMaskLowercase, const uchar letterCode[],
		 const std::string &mainSequenceAlphabet);

  void toFiles( const std::string& baseName,
		bool isAppendPrj, size_t textLength ) const;

  // Find the smallest match to the text, starting at the given
  // position in the query, such that there are at most maxHits
  // matches, and the match-depth is at least minDepth, or the
  // match-depth is maxDepth.  Return the range of matching indices
  // via begPtr and endPtr.
  void match( const PosPart *&begPtr, const PosPart *&endPtr,
              const uchar* queryPtr, const uchar* text, unsigned seedNum,
              size_t maxHits, size_t minDepth, size_t maxDepth ) const;

  // Count matches of all sizes (up to maxDepth), starting at the
  // given position in the query.
  void countMatches( std::vector<unsigned long long>& counts,
		     const uchar* queryPtr, const uchar* text,
		     unsigned seedNum, size_t maxDepth ) const;

private:
  std::vector<CyclicSubsetSeed> seeds;
  std::vector<const OffPart *> bucketEnds;
  std::vector<const indexT *> bucketStepEnds;

  VectorOrMmap<PosPart> suffixArray;  // sorted indices
  VectorOrMmap<OffPart> buckets;
  std::vector<indexT> bucketSteps;  // step size for each k-mer

  VectorOrMmap<indexT> childTable;
  VectorOrMmap<unsigned short> kiddyTable;  // smaller child table
  VectorOrMmap<unsigned char> chibiTable;  // even smaller child table

  enum ChildDirection { FORWARD, REVERSE, UNKNOWN };

  // This does the same thing as equalRange, but uses a child table:
  void childRange( indexT& beg, indexT& end, ChildDirection& childDirection,
                   const uchar* textBase,
                   const uchar* subsetMap, uchar subset ) const;

  // Return the maximum prefix size covered by the buckets.
  size_t maxBucketPrefix(unsigned seedNum) const
  { return bucketStepEnds[seedNum + 1] - bucketStepEnds[seedNum] - 1; }

  void makeBucketSteps(const unsigned *bucketDepths, size_t wordLength);

  size_t bucketsSize() const {
    size_t n = offParts;
    for (size_t i = 0; i < seeds.size(); ++i) {
      n += bucketStepEnds[i][0];
    }
    return n;
  }

  void initBucketEnds() {
    bucketEnds.resize(seeds.size());
    const OffPart *p = &buckets[0];
    for (size_t i = 0; i < seeds.size(); ++i) {
      bucketEnds[i] = p;
      p += bucketStepEnds[i][0];
    }
  }

  void sort2( const uchar* text, const CyclicSubsetSeed& seed,
	      PosPart *beg, const uchar* subsetMap );

  void radixSort1( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth );
  void radixSort2( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth );
  void radixSort3( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth );
  void radixSort4( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth );
  void radixSortN( std::vector<Range>& rangeStack,
		   const uchar* text, const uchar* subsetMap,
		   PosPart *beg, PosPart *end, indexT depth,
		   unsigned subsetCount, indexT* bucketSize );

  void sortRanges( std::vector<Range>* stacks, indexT* bucketSizes,
		   const uchar* text,
		   unsigned wordLength, const CyclicSubsetSeed& seed,
		   size_t maxUnsortedInterval, size_t numOfThreads );

  indexT getChildForward( indexT from ) const{
    return
      !childTable.empty() ? childTable[ from ] :
      !kiddyTable.empty() ? from + kiddyTable[ from ] :
      !chibiTable.empty() ? from + chibiTable[ from ] : from;
  }

  indexT getChildReverse( indexT from ) const{
    return
      !childTable.empty() ? childTable[ from - 1 ] :
      !kiddyTable.empty() ? from - kiddyTable[ from - 1 ] :
      !chibiTable.empty() ? from - chibiTable[ from - 1 ] : from;
  }

  void setKiddy( indexT index, indexT value ){
    kiddyTable.v[ index ] = (value < USHRT_MAX) ? value : 0;
  }

  void setChibi( indexT index, indexT value ){
    chibiTable.v[ index ] = (value < UCHAR_MAX) ? value : 0;
  }

  void setChildForward(const PosPart *from, const PosPart *to) {
    if( to == from ) return;
    const PosPart *origin = &suffixArray.v[0];
    indexT i = posCount(origin, from);
    /**/ if (!childTable.v.empty()) childTable.v[i] = posCount(origin, to);
    else if (!kiddyTable.v.empty()) setKiddy(i, posCount(from, to));
    else if (!chibiTable.v.empty()) setChibi(i, posCount(from, to));
  }

  void setChildReverse(const PosPart *from, const PosPart *to) {
    if( to == from ) return;
    const PosPart *origin = &suffixArray.v[0];
    indexT i = posCount(origin, from) - 1;
    /**/ if (!childTable.v.empty()) childTable.v[i] = posCount(origin, to);
    else if (!kiddyTable.v.empty()) setKiddy(i, posCount(to, from));
    else if (!chibiTable.v.empty()) setChibi(i, posCount(to, from));
  }

  bool isChildDirectionForward(const PosPart *beg) const {
    const PosPart *origin = &suffixArray.v[0];
    indexT i = posCount(origin, beg);
    return
      !childTable.v.empty() ? childTable.v[ i ] == 0 :
      !kiddyTable.v.empty() ? kiddyTable.v[ i ] == USHRT_MAX :
      !chibiTable.v.empty() ? chibiTable.v[ i ] == UCHAR_MAX : true;
  }
};

}  // end namespace
#endif
