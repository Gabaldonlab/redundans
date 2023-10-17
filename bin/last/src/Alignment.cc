// Copyright 2008, 2009, 2011, 2012, 2013, 2014 Martin C. Frith

#include "Alignment.hh"
#include "Alphabet.hh"
#include "GeneticCode.hh"
#include "TwoQualityScoreMatrix.hh"

#include <assert.h>

// make C++ tolerable:
#define IT(type) std::vector<type>::iterator

using namespace cbrc;

static void addSeedCounts(const uchar *seq1, const uchar *seq2, size_t size,
			  double *counts) {
  for (size_t i = 0; i < size; ++i) {
    ++counts[seq1[i] * scoreMatrixRowSize + seq2[i]];
  }
  counts[scoreMatrixRowSize * scoreMatrixRowSize] += size;
}

// Does x precede and touch y in both sequences?
static bool isNext( const SegmentPair& x, const SegmentPair& y ){
  return x.end1() == y.beg1() && x.end2() == y.beg2();
}

void Alignment::makeXdrop( Aligners &aligners, bool isGreedy, bool isFullScore,
			   const uchar* seq1, const uchar* seq2, int globality,
			   const ScoreMatrixRow* scoreMatrix,
			   int smMax, int smMin,
			   const const_dbl_ptr* probMatrix, double scale,
			   const GapCosts& gap, int maxDrop,
			   size_t frameSize, const ScoreMatrixRow* pssm2,
                           const TwoQualityScoreMatrix& sm2qual,
                           const uchar* qual1, const uchar* qual2,
			   const Alphabet& alph, AlignmentExtras& extras,
			   double gamma, int outputType ){
  if (probMatrix) score = seed.score;  // else keep the old score
  if (outputType > 3 && !isFullScore) extras.fullScore = seed.score;

  if( outputType == 7 ){
    assert( seed.size > 0 );  // makes things easier to understand
    const int numOfTransitions = frameSize ? 9 : 5;
    std::vector<double> &ec = extras.expectedCounts;
    ec.resize(scoreMatrixRowSize * scoreMatrixRowSize + numOfTransitions);
    addSeedCounts(seq1 + seed.beg1(), seq2 + seed.beg2(), seed.size, &ec[0]);
  }

  // extend a gapped alignment in the left/reverse direction from the seed:
  blocks.clear();
  std::vector<char>& columnAmbiguityCodes = extras.columnAmbiguityCodes;
  extend( blocks, columnAmbiguityCodes, aligners, isGreedy, isFullScore,
	  seq1, seq2, seed.beg1(), seed.beg2(), false, globality,
	  scoreMatrix, smMax, smMin, probMatrix, scale, maxDrop, gap,
	  frameSize, pssm2, sm2qual, qual1, qual2, alph,
	  extras, gamma, outputType );

  if( score == -INF ) return;  // maybe unnecessary?

  // convert left-extension coordinates to sequence coordinates:
  SegmentPair::indexT seedBeg1 = seed.beg1();
  SegmentPair::indexT seedBeg2 = aaToDna( seed.beg2(), frameSize );
  for( IT(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    i->start1 = seedBeg1 - i->start1 - i->size;
    // careful: i->start2 might be -1 (reverse frameshift)
    i->start2 = dnaToAa( seedBeg2 - i->start2, frameSize ) - i->size;
  }

  // extend a gapped alignment in the right/forward direction from the seed:
  std::vector<SegmentPair> forwardBlocks;
  std::vector<char> forwardAmbiguities;
  extend( forwardBlocks, forwardAmbiguities, aligners, isGreedy, isFullScore,
	  seq1, seq2, seed.end1(), seed.end2(), true, globality,
	  scoreMatrix, smMax, smMin, probMatrix, scale, maxDrop, gap,
	  frameSize, pssm2, sm2qual, qual1, qual2, alph,
	  extras, gamma, outputType );

  if( score == -INF ) return;  // maybe unnecessary?

  // convert right-extension coordinates to sequence coordinates:
  SegmentPair::indexT seedEnd1 = seed.end1();
  SegmentPair::indexT seedEnd2 = aaToDna( seed.end2(), frameSize );
  for( IT(SegmentPair) i = forwardBlocks.begin(); i < forwardBlocks.end();
       ++i ){
    i->start1 = seedEnd1 + i->start1;
    // careful: i->start2 might be -1 (reverse frameshift)
    i->start2 = dnaToAa( seedEnd2 + i->start2, frameSize );
  }

  bool isMergeSeedReverse = !blocks.empty() && isNext( blocks.back(), seed );
  bool isMergeSeedForward =
    !forwardBlocks.empty() && isNext( seed, forwardBlocks.back() );

  if( seed.size == 0 && !isMergeSeedReverse && !isMergeSeedForward ){
    // unusual, weird case: give up
    score = -INF;
    return;
  }

  // splice together the two extensions and the seed (a bit messy):

  blocks.reserve( blocks.size() + forwardBlocks.size() +
		  1 - isMergeSeedReverse - isMergeSeedForward );

  if( isMergeSeedReverse ) blocks.back().size += seed.size;
  else                     blocks.push_back(seed);

  if( isMergeSeedForward ){
    blocks.back().size += forwardBlocks.back().size;
    forwardBlocks.pop_back();
  }

  size_t oldSize = blocks.size();
  blocks.insert( blocks.end(), forwardBlocks.rbegin(), forwardBlocks.rend() );
  for (size_t i = oldSize; i < blocks.size(); ++i)
    blocks[i - 1].score = blocks[i].score;

  if( outputType > 3 ){  // set the un-ambiguity of the core to a max value:
    columnAmbiguityCodes.insert( columnAmbiguityCodes.end(), seed.size, 126 );
  }

  columnAmbiguityCodes.insert( columnAmbiguityCodes.end(),
                               forwardAmbiguities.rbegin(),
                               forwardAmbiguities.rend() );
}

// cost of the gap between x and y
static int gapCost(const SegmentPair &x, const SegmentPair &y,
		   const GapCosts &gapCosts, size_t frameSize) {
  if (gapCosts.isNewFrameshifts()) return x.score;
  size_t gapSize1 = y.beg1() - x.end1();
  size_t gapSize2, frameshift2;
  sizeAndFrameshift(x.end2(), y.beg2(), frameSize, gapSize2, frameshift2);
  int cost = gapCosts.cost(gapSize1, gapSize2);
  if (frameshift2) cost += gapCosts.frameshiftCost;
  return cost;
}

bool Alignment::isOptimal( const uchar* seq1, const uchar* seq2, int globality,
			   const ScoreMatrixRow* scoreMatrix, int maxDrop,
			   const GapCosts& gapCosts, size_t frameSize,
			   const ScoreMatrixRow* pssm2,
                           const TwoQualityScoreMatrix& sm2qual,
                           const uchar* qual1, const uchar* qual2 ) const{
  int maxScore = 0;
  int runningScore = 0;

  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];

    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];
      runningScore -= gapCost( x, y, gapCosts, frameSize );
      if( !globality && runningScore <= 0 ) return false;
      if( runningScore < maxScore - maxDrop ) return false;
    }

    const uchar* s1 = seq1 + y.beg1();
    const uchar* s2 = seq2 + y.beg2();
    const uchar* e1 = seq1 + y.end1();

    const ScoreMatrixRow* p2 = pssm2 ? pssm2 + y.beg2() : 0;
    const uchar* q1 = qual1 ? qual1 + y.beg1() : 0;
    const uchar* q2 = qual2 ? qual2 + y.beg2() : 0;

    while( s1 < e1 ){
      /**/ if( sm2qual ) runningScore += sm2qual( *s1++, *s2++, *q1++, *q2++ );
      else if( pssm2 )   runningScore += ( *p2++ )[ *s1++ ];
      else               runningScore += scoreMatrix[ *s1++ ][ *s2++ ];

      if( runningScore > maxScore ) maxScore = runningScore;
      else if( !globality && runningScore <= 0 ) return false;
      else if( !globality && (s1 == e1 && i+1 == blocks.size()) ) return false;
      else if( runningScore < maxScore - maxDrop ) return false;
    }
  }

  return true;
}

bool Alignment::hasGoodSegment(const uchar *seq1, const uchar *seq2,
			       int minScore, const ScoreMatrixRow *scoreMatrix,
			       const GapCosts &gapCosts, size_t frameSize,
			       const ScoreMatrixRow *pssm2,
			       const TwoQualityScoreMatrix &sm2qual,
			       const uchar *qual1, const uchar *qual2) const {
  int score = 0;

  for (size_t i = 0; i < blocks.size(); ++i) {
    const SegmentPair& y = blocks[i];

    if (i > 0) {  // between each pair of aligned blocks:
      score -= gapCost(blocks[i - 1], y, gapCosts, frameSize);
      if (score < 0) score = 0;
    }

    const uchar *s1 = seq1 + y.beg1();
    const uchar *s2 = seq2 + y.beg2();
    const uchar *e1 = seq1 + y.end1();

    const ScoreMatrixRow *p2 = pssm2 ? pssm2 + y.beg2() : 0;
    const uchar *q1 = qual1 ? qual1 + y.beg1() : 0;
    const uchar *q2 = qual2 ? qual2 + y.beg2() : 0;

    while (s1 < e1) {
      /**/ if (sm2qual) score += sm2qual(*s1++, *s2++, *q1++, *q2++);
      else if (pssm2)   score += (*p2++)[*s1++];
      else              score += scoreMatrix[*s1++][*s2++];

      if (score >= minScore) return true;
      if (score < 0) score = 0;
    }
  }

  return false;
}

static void getColumnCodes(const Centroid& centroid, std::vector<char>& codes,
			   const std::vector<SegmentPair>& chunks,
			   bool isForward) {
  for (size_t i = 0; i < chunks.size(); ++i) {
    const SegmentPair& x = chunks[i];
    centroid.getMatchAmbiguities(codes, x.end1(), x.end2(), x.size);
    size_t j = i + 1;
    bool isNext = (j < chunks.size());
    size_t end1 = isNext ? chunks[j].end1() : 0;
    size_t end2 = isNext ? chunks[j].end2() : 0;
    // ASSUMPTION: if there is an insertion adjacent to a deletion,
    // the deletion will get printed first.
    if (isForward) {
      centroid.getInsertAmbiguities(codes, x.beg2(), end2);
      centroid.getDeleteAmbiguities(codes, x.beg1(), end1);
    } else {
      centroid.getDeleteAmbiguities(codes, x.beg1(), end1);
      centroid.getInsertAmbiguities(codes, x.beg2(), end2);
    }
  }
}

static void getColumnCodes(const FrameshiftXdropAligner &fxa,
			   std::vector<char> &codes,
			   const std::vector<SegmentPair> &chunks) {
  for (size_t i = 0; i < chunks.size(); ++i) {
    const SegmentPair &x = chunks[i];
    for (size_t k = x.size; k-- > 0;) {
      double p = fxa.matchProb(x.beg1() + k, x.beg2() + k * 3);
      codes.push_back(asciiProbability(p));
    }
    size_t j = i + 1;
    bool isNext = (j < chunks.size());
    size_t end1 = isNext ? chunks[j].end1() : 0;
    size_t end2 = isNext ? chunks[j].beg2() + chunks[j].size * 3 : 0;
    size_t n1 = x.beg1() - end1;
    size_t n2 = (x.beg2() - end2 + 1) / 3;
    codes.insert(codes.end(), n1 + n2, '-');
  }
}

void Alignment::extend( std::vector< SegmentPair >& chunks,
			std::vector< char >& columnCodes,
			Aligners &aligners, bool isGreedy, bool isFullScore,
			const uchar* seq1, const uchar* seq2,
			size_t start1, size_t start2,
			bool isForward, int globality,
			const ScoreMatrixRow* sm, int smMax, int smMin,
			const const_dbl_ptr* probMat, double scale,
			int maxDrop, const GapCosts& gap, size_t frameSize,
			const ScoreMatrixRow* pssm2,
			const TwoQualityScoreMatrix& sm2qual,
                        const uchar* qual1, const uchar* qual2,
			const Alphabet& alph, AlignmentExtras& extras,
			double gamma, int outputType ){
  const GapCosts::Piece &del = gap.delPieces[0];
  const GapCosts::Piece &ins = gap.insPieces[0];
  Centroid &centroid = aligners.centroid;
  GappedXdropAligner& aligner = centroid.aligner();
  GreedyXdropAligner &greedyAligner = aligners.greedyAligner;

  double *subsCounts[scoreMatrixRowSize];
  double *tranCounts;
  if (outputType == 7) {
    double *ec = &extras.expectedCounts[0];
    for (int i = 0; i < scoreMatrixRowSize; ++i)
      subsCounts[i] = ec + i * scoreMatrixRowSize;
    tranCounts = ec + scoreMatrixRowSize * scoreMatrixRowSize;
  }

  if( frameSize ){
    assert( !isGreedy );
    assert( !globality );
    assert( !pssm2 );
    assert( !sm2qual );

    size_t dnaStart = aaToDna( start2, frameSize );
    size_t frame1 = isForward ? dnaStart + 1 : dnaStart - 1;
    size_t end1, end2, size;
    int gapCost;

    if (gap.isNewFrameshifts()) {
      assert(isFullScore);
      size_t frame2 = isForward ? dnaStart + 2 : dnaStart - 2;
      aligner.alignFrame(seq1 + start1, seq2 + start2,
			 seq2 + dnaToAa(frame1, frameSize),
			 seq2 + dnaToAa(frame2, frameSize),
			 isForward, sm, gap, maxDrop);
      while (aligner.getNextChunkFrame(end1, end2, size, gapCost, gap))
	chunks.push_back(SegmentPair(end1 - size, end2 - size * 3, size,
				     gapCost));
      if (!probMat) return;
      FrameshiftXdropAligner &fxa = aligners.frameshiftAligner;
      double probDropLimit = exp(scale * -maxDrop);
      double s = fxa.forward(seq1 + start1, seq2 + start2,
			     seq2 + dnaToAa(frame1, frameSize),
			     seq2 + dnaToAa(frame2, frameSize),
			     isForward, probMat, gap, probDropLimit);
      score += s / scale;
      if (outputType < 4) return;
      fxa.backward(isForward, probMat, gap);
      getColumnCodes(fxa, columnCodes, chunks);
      if (outputType == 7) fxa.count(isForward, gap, subsCounts, tranCounts);
    } else {
      assert(!isFullScore);
      assert(outputType < 4);
      size_t frame2 = isForward ? dnaStart - 1 : dnaStart + 1;
      score += aligner.align3(seq1 + start1, seq2 + start2,
			      seq2 + dnaToAa(frame1, frameSize),
			      seq2 + dnaToAa(frame2, frameSize), isForward,
			      sm, del.openCost, del.growCost, gap.pairCost,
			      gap.frameshiftCost, maxDrop, smMax);
      // This should be OK even if end2 < size * 3:
      while (aligner.getNextChunk3(end1, end2, size,
				   del.openCost, del.growCost, gap.pairCost,
				   gap.frameshiftCost))
	chunks.push_back(SegmentPair(end1 - size, end2 - size * 3, size));
    }

    return;
  }

  if (!isForward) {
    --start1;
    --start2;
  }
  seq1 += start1;
  seq2 += start2;

  bool isSimdMatrix = (alph.size == 4 && !globality && gap.isAffine &&
		       smMin >= SCHAR_MIN &&
		       maxDrop + smMax * 2 - smMin < UCHAR_MAX);
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      if (sm[i][j] != sm[alph.numbersToLowercase[i]][j])
	isSimdMatrix = false;

  int extensionScore =
    isGreedy  ? greedyAligner.align(seq1, seq2, isForward, sm,
				    maxDrop, alph.size)
    : sm2qual ? aligner.align2qual(seq1, qual1 + start1, seq2, qual2 + start2,
				   isForward, globality, sm2qual,
				   del.openCost, del.growCost,
				   ins.openCost, ins.growCost,
				   gap.pairCost, gap.isAffine, maxDrop, smMax)
    : pssm2   ? aligner.alignPssm(seq1, pssm2 + start2, isForward, globality,
				  del.openCost, del.growCost,
				  ins.openCost, ins.growCost,
				  gap.pairCost, gap.isAffine, maxDrop, smMax)
#if defined __SSE4_1__ || defined __ARM_NEON
    : isSimdMatrix ? aligner.alignDna(seq1, seq2, isForward, sm,
				      del.openCost, del.growCost,
				      ins.openCost, ins.growCost,
				      maxDrop, smMax, alph.numbersToUppercase)
#endif
    :           aligner.align(seq1, seq2, isForward, globality, sm,
			      del.openCost, del.growCost,
			      ins.openCost, ins.growCost,
			      gap.pairCost, gap.isAffine, maxDrop, smMax);

  if( extensionScore == -INF ){
    score = -INF;  // avoid score overflow
    return;  // avoid ill-defined probabilistic alignment
  }

  if( outputType < 5 || outputType > 6 ){  // ordinary max-score alignment
    size_t end1, end2, size;
    if( isGreedy ){
      while( greedyAligner.getNextChunk( end1, end2, size ) )
	chunks.push_back( SegmentPair( end1 - size, end2 - size, size ) );
    }
#if defined __SSE4_1__ || defined __ARM_NEON
    else if (isSimdMatrix && !pssm2 && !sm2qual) {
      while (aligner.getNextChunkDna(end1, end2, size,
				     del.openCost, del.growCost,
				     ins.openCost, ins.growCost))
	chunks.push_back(SegmentPair(end1 - size, end2 - size, size));
    }
#endif
    else {
      while( aligner.getNextChunk( end1, end2, size,
				   del.openCost, del.growCost,
				   ins.openCost, ins.growCost, gap.pairCost ) )
	chunks.push_back( SegmentPair( end1 - size, end2 - size, size ) );
    }
  }

  if (!probMat) return;
  if (!isFullScore) score += extensionScore;

  if (outputType > 3 || isFullScore) {
    assert( !isGreedy );
    assert( !sm2qual );
    double s = centroid.forward(seq1, seq2, start2, isForward,
				probMat, gap, globality);
    if (isFullScore) {
      score += s / scale;
    } else {
      extras.fullScore += s / scale;
    }
    if (outputType < 4) return;
    centroid.backward(isForward, probMat, gap, globality);
    if (outputType > 4 && outputType < 7) {  // gamma-centroid / LAMA alignment
      centroid.dp(outputType, gamma);
      centroid.traceback(chunks, outputType, gamma);
    }
    getColumnCodes(centroid, columnCodes, chunks, isForward);
    if (outputType == 7) {
      centroid.addExpectedCounts(start2, isForward, probMat, gap, alph.size,
				 subsCounts, tranCounts);
    }
  }
}
