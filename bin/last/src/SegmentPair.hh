// Copyright 2008, 2009, 2010, 2011 Martin C. Frith

// This struct holds a pair of equal-length segments, in a pair of
// sequences.  In other words, it holds a gapless alignment.

#ifndef SEGMENT_PAIR_HH
#define SEGMENT_PAIR_HH

#include <stddef.h>

namespace cbrc{

struct SegmentPair{
  typedef unsigned char uchar;

#if LAST_POS_BYTES > 4
  typedef size_t indexT;
#else
  typedef unsigned indexT;
#endif

  SegmentPair(){}

  SegmentPair( indexT s1, indexT s2, indexT sz, int sc = 0 )
    : start1(s1), start2(s2), size(sz), score(sc){}

  // Shrink the SegmentPair to the longest run of identical letters
  // within it.  Define "identical" as:
  // map1[ letter from seq1 ] == map2[ letter from seq2 ].
  void maxIdenticalRun( const uchar* seq1, const uchar* seq2,
			const uchar* map1, const uchar* map2 );

  indexT beg1() const{ return start1; }         // start in sequence 1
  indexT beg2() const{ return start2; }         // start in sequence 2
  indexT end1() const{ return start1 + size; }  // end in sequence 1
  indexT end2() const{ return start2 + size; }  // end in sequence 2
  indexT diagonal() const{ return start1 - start2; }  // may wrap around!

  bool operator==( const SegmentPair& s ) const{
    return start1 == s.start1
        && start2 == s.start2
        && size == s.size
        && score == s.score; }

  indexT start1;
  indexT start2;
  indexT size;
  int score;
};

}

#endif
