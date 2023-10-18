// Copyright 2016 Martin C. Frith

#include "SubsetMinimizerFinder.hh"
#include "CyclicSubsetSeed.hh"

namespace cbrc {

void SubsetMinimizerFinder::init(const CyclicSubsetSeed &seed,
				 const uchar *beg,
				 const uchar *end) {
  const uchar *m = seed.firstMap();
  while (beg < end && m[*beg] == CyclicSubsetSeed::DELIMITER) ++beg;
  minima.assign(1, beg);
}

bool SubsetMinimizerFinder::isMinimizer(const CyclicSubsetSeed &seed,
					const uchar *pos,
					const uchar *end,
					size_t window) {
  const uchar *subsetMap = seed.firstMap();

  while (true) {
    const uchar *currentMinimum = minima[0];
    if (currentMinimum > pos) return false;
    if (currentMinimum == pos) return true;
    const uchar *newPos = minima.back() + 1;
    if (newPos == end) return false;
    if (subsetMap[*newPos] == CyclicSubsetSeed::DELIMITER) {
      init(seed, newPos + 1, end);
      continue;
    }
    size_t diff = newPos - currentMinimum;
    size_t stop = (diff >= window);
    size_t i = minima.size();
    while (i > stop && seed.isLess(newPos, minima[i - 1], subsetMap)) {
      --i;
    }
    minima.resize(i);
    if (stop) minima.erase(minima.begin());
    minima.push_back(newPos);
  }
}

}
