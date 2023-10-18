// Copyright 2008, 2010 Martin C. Frith

#include "AlignmentPot.hh"

namespace cbrc{

void AlignmentPot::eraseSuboptimal(){
  if( items.empty() ) return;

  std::sort( items.begin(), items.end(), lessBeg );

  for( iterator i = items.begin() + 1; i < items.end(); ++i ){
    if( i->beg1() == (i-1)->beg1() && i->beg2() == (i-1)->beg2() ){
      mark( *i );
    }
  }

  std::sort( items.begin(), items.end(), lessEnd );

  for( iterator i = items.begin() + 1; i < items.end(); ++i ){
    if( i->end1() == (i-1)->end1() && i->end2() == (i-1)->end2() ){
      mark( *i );
    }
  }

  items.erase( std::remove_if( items.begin(), items.end(), isMarked ),
	       items.end() );
}

}  // end namespace cbrc
