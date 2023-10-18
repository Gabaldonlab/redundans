// Copyright 2008 Martin C. Frith

// This struct holds alignments, and includes a procedure to
// non-redundantize alignments that share endpoints.

#ifndef ALIGNMENTPOT_HH
#define ALIGNMENTPOT_HH
#include "Alignment.hh"
#include <algorithm>  // sort
#include <vector>

namespace cbrc{

struct AlignmentPot{
  typedef std::vector<Alignment>::iterator iterator;

  // add an alignment to the pot
  void add( const Alignment& aln ) { items.push_back(aln); }

  // the number of alignments in the pot
  size_t size() const { return items.size(); }

  // if several alignments share an endpoint, erase all but one
  // highest-scoring alignment
  void eraseSuboptimal();

  // sort the alignments in descending order of score
  void sort() { std::sort( items.begin(), items.end(), moreScore ); }

  // data:
  std::vector<Alignment> items;

  static bool moreScore( const Alignment& x, const Alignment& y ){
    // Try to break ties, so that alignments come in a consistent
    // order.  This makes it easier to compare different results.
    return x.score != y.score ? x.score > y.score : lessBeg( x, y );
  }

  static bool lessBeg( const Alignment& x, const Alignment& y ){
    if( x.beg1() != y.beg1() ) return x.beg1() < y.beg1();
    if( x.beg2() != y.beg2() ) return x.beg2() < y.beg2();
    if( x.score  != y.score  ) return x.score  > y.score;
    // arbitrary (but systematic) order:
    const SegmentPair& a = x.seed;
    const SegmentPair& b = y.seed;
    return a.beg1() != b.beg1() ? a.beg1() < b.beg1() : a.beg2() < b.beg2();
  }

  static bool lessEnd( const Alignment& x, const Alignment& y ){
    if( x.end1() != y.end1() ) return x.end1() < y.end1();
    if( x.end2() != y.end2() ) return x.end2() < y.end2();
    if( x.score  != y.score  ) return x.score  > y.score;
    // arbitrary (but systematic) order:
    const SegmentPair& a = x.seed;
    const SegmentPair& b = y.seed;
    return a.beg1() != b.beg1() ? a.beg1() < b.beg1() : a.beg2() < b.beg2();
  }

  static void mark( Alignment& a ) { a.blocks[0].score = -1; }

  static bool isMarked( Alignment& a ) { return a.blocks[0].score == -1; }
};

}  // end namespace cbrc
#endif  // ALIGNMENTPOT_HH
