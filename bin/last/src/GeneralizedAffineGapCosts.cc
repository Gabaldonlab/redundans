// Copyright 2008, 2009, 2012 Martin C. Frith

#include "GeneralizedAffineGapCosts.hh"
#include <algorithm>

namespace cbrc{

int GeneralizedAffineGapCosts::cost( int gapSize1, int gapSize2 ) const{
  int delPart = gapSize1 ? delExist + delExtend * gapSize1 : 0;
  int insPart = gapSize2 ? insExist + insExtend * gapSize2 : 0;
  int c = delPart + insPart;

  if( gapSize1 >= gapSize2 && pairExtend < insExist + insExtend + delExtend ){
    int d = delPart + (pairExtend - delExtend) * gapSize2;
    c = std::min( c, d );
  }

  if( gapSize2 >= gapSize1 && pairExtend < delExist + delExtend + insExtend ){
    int d = insPart + (pairExtend - insExtend) * gapSize1;
    c = std::min( c, d );
  }

  return c;
}

}
