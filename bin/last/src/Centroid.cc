// Copyright 2008, 2009, 2010, 2011 Michiaki Hamada
// Copyright 2012, 2013 Toshiyuki Sato

#include "Centroid.hh"
#include "GappedXdropAlignerInl.hh"

#include <algorithm>
#include <cmath> // for exp
#include <cfloat>   // for DBL_MAX

static const double DINF = DBL_MAX / 2;

namespace{
  double EXP ( double x ) {
    return std::exp (x);
  }
}

namespace cbrc{

  void Centroid::setPssm( const ScoreMatrixRow* pssm, size_t qsize, double T,
			  const OneQualityExpMatrix& oqem,
			  const uchar* sequenceBeg, const uchar* qualityBeg ) {
    pssmExp.resize( qsize * scoreMatrixRowSize );
    pssmExp2 = reinterpret_cast<ExpMatrixRow*> ( &pssmExp[0] );

    if( oqem ){  // fast special case
      makePositionSpecificExpMatrix( oqem, sequenceBeg, sequenceBeg + qsize,
                                     qualityBeg, &pssmExp[0] );
    }
    else{  // slow general case
      for ( size_t i=0; i<qsize; ++i ) {
        for ( unsigned j=0; j<scoreMatrixRowSize; ++j ) {
          pssmExp2[ i ][ j ] = EXP ( pssm[ i ][ j ] / T );
        }
      }
    }
  }

  void Centroid::setLetterProbsPerPosition(unsigned alphabetSize,
					   size_t sequenceLength,
					   const uchar *sequence,
					   const uchar *qualityCodes,
					   bool isFastq,
					   const double *qualToProbCorrect,
					   const double *letterProbs,
					   const uchar *toUnmasked) {
    letterProbsPerPosition.resize(sequenceLength * alphabetSize);
    for (size_t i = 0; i < sequenceLength; ++i) {
      size_t j = i * alphabetSize;
      double *out = &letterProbsPerPosition[j];
      if (isFastq) {
	unsigned letter = toUnmasked[sequence[i]];
	if (letter < alphabetSize) {
	  double p = qualToProbCorrect[qualityCodes[i]];
	  for (unsigned k = 0; k < alphabetSize; ++k) {
	    out[k] = (1 - p) * letterProbs[k];
	  }
	  out[letter] = p * (1 - letterProbs[letter]);
	  // it's OK to scale the "out" values by a constant, per "i"
	} else {
	  std::fill_n(out, alphabetSize, 0.0);
	}
      } else {
	for (unsigned k = 0; k < alphabetSize; ++k) {
	  out[k] = qualToProbCorrect[qualityCodes[j + k]];
	}
      }
    }
  }

  double Centroid::forward(const uchar *seq1, const uchar *seq2,
			   size_t start2, bool isExtendFwd,
			   const const_dbl_ptr *substitutionProbs,
			   const GapCosts &gapCosts, int globality) {
    seq1ptr = seq1;
    seq2ptr = seq2;
    pssmPtr = pssmExp.empty() ? 0 : pssmExp2 + start2;
    const int seqIncrement = isExtendFwd ? 1 : -1;

    const double delInit = gapCosts.delProbPieces[0].openProb;
    const double delNext = gapCosts.delProbPieces[0].growProb;
    const double insInit = gapCosts.insProbPieces[0].openProb;
    const double insNext = gapCosts.insProbPieces[0].growProb;

    initForward();

    size_t antidiagonal = 0;
    size_t seq1beg = 0;
    size_t diagPos = xdropPadLen - 1;
    size_t horiPos = xdropPadLen * 2 - 1;
    size_t thisPos = xdropPadLen * 2;

    double sumOfEdgeProbRatios = 0;
    double sumOfProbRatios = 0;
    double logSumOfProbRatios = 0;

    while (1) {
      double *fM0 = &fM[thisPos];
      double *fD0 = &fD[thisPos];
      double *fI0 = &fI[thisPos];
      for (int i = 0; i < xdropPadLen; ++i) {
	*fM0++ = *fD0++ = *fI0++ = 0.0;
      }
      thisPos += xdropPadLen;

      const double *fD1 = &fD[horiPos];
      const double *fI1 = &fI[horiPos + 1];
      const double *fM2 = &fM[diagPos];

      ++antidiagonal;
      const size_t nextPos = xa.scoreEndIndex(antidiagonal);
      const int numCells = nextPos - thisPos;
      const uchar *s1 = seq1ptr;

      if (!pssmPtr) {
	const uchar *s2 = seq2ptr;
	for (int i = 0; i < numCells; ++i) {
	  const double matchProb = substitutionProbs[*s1][*s2];
	  const double xD = fD1[i];
	  const double xI = fI1[i];
	  const double xSum = fM2[i] + xD + xI;
	  fD0[i] = xSum * delInit + xD * delNext;
	  fI0[i] = xSum * insInit + xI * insNext;
	  fM0[i] = xSum * matchProb;
	  sumOfProbRatios += xSum;
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}
      } else {
	const ExpMatrixRow *p2 = pssmPtr;
	for (int i = 0; i < numCells; ++i) {
	  const double matchProb = (*p2)[*s1];
	  const double xD = fD1[i];
	  const double xI = fI1[i];
	  const double xSum = fM2[i] + xD + xI;
	  fD0[i] = xSum * delInit + xD * delNext;
	  fI0[i] = xSum * insInit + xI * insNext;
	  fM0[i] = xSum * matchProb;
	  sumOfProbRatios += xSum;
	  s1 += seqIncrement;
	  p2 -= seqIncrement;
	}
      }

      if (globality) {
	const int n = numCells - 1;
	if (substitutionProbs[0][*seq2ptr] <= 0) {
	  sumOfEdgeProbRatios += fM2[0] + fD1[0] + fI1[0];
	}
	if (n > 0 && substitutionProbs[*(s1 - seqIncrement)][0] <= 0) {
	  sumOfEdgeProbRatios += fM2[n] + fD1[n] + fI1[n];
	}
      }

      if (antidiagonal == numAntidiagonals) break;

      diagPos = horiPos;
      horiPos = thisPos - 1;
      thisPos = nextPos;

      const size_t newSeq1beg = xa.seq1start(antidiagonal);
      if (newSeq1beg > seq1beg) {
	seq1beg = newSeq1beg;
	seq1ptr += seqIncrement;
	++diagPos;
	++horiPos;
      } else {
	seq2ptr += seqIncrement;
	if (pssmPtr) pssmPtr += seqIncrement;
      }

      if (antidiagonal % rescaleStep == 0) {
	const double scale = 1 / sumOfProbRatios;
	rescales[antidiagonal / rescaleStep - 1] = scale;
	rescaleFwdProbs(xa.scoreEndIndex(antidiagonal - 2), thisPos, scale);
	logSumOfProbRatios += log(sumOfProbRatios);
	sumOfEdgeProbRatios *= scale;
	sumOfProbRatios = 1;
      }
    }

    if (globality) {
      assert(sumOfEdgeProbRatios > 0);
      sumOfProbRatios = sumOfEdgeProbRatios;
    }
    rescaledSumOfProbRatios = sumOfProbRatios;
    return logSumOfProbRatios + log(sumOfProbRatios);
  }

  // added by M. Hamada
  // compute posterior probabilities while executing backward algorithm
  void Centroid::backward(bool isExtendFwd,
			  const const_dbl_ptr *substitutionProbs,
			  const GapCosts& gapCosts, int globality) {
    const int seqIncrement = isExtendFwd ? 1 : -1;

    const double delInit = gapCosts.delProbPieces[0].openProb;
    const double delNext = gapCosts.delProbPieces[0].growProb;
    const double insInit = gapCosts.insProbPieces[0].openProb;
    const double insNext = gapCosts.insProbPieces[0].growProb;

    size_t antidiagonal = numAntidiagonals - 1;
    size_t seq1beg = xa.seq1start(antidiagonal);
    size_t oldPos = xa.scoreEndIndex(numAntidiagonals);
    initBackward(oldPos);
    double scaledUnit = 1 / rescaledSumOfProbRatios;

    while (1) {
      const size_t newPos = xa.scoreEndIndex(antidiagonal);
      const double *bM0 = &bM[newPos + xdropPadLen];
      const double *bD0 = &bD[newPos + xdropPadLen];
      const double *bI0 = &bI[newPos + xdropPadLen];

      const double *fD0 = &fD[newPos + xdropPadLen];
      const double *fI0 = &fI[newPos + xdropPadLen];

      const size_t vertPos = xa.vert(antidiagonal, seq1beg);
      const size_t diagPos = xa.diag(antidiagonal, seq1beg);
      double *bD1 = &bD[vertPos - 1];
      double *bI1 = &bI[vertPos];
      double *bM2 = &bM[diagPos];

      const size_t seq2pos = antidiagonal - seq1beg;
      double *mDout = &mD[seq1beg + 1];
      double *mIout = &mI[seq2pos + 1];

      const int numCells = oldPos - newPos - xdropPadLen;
      const uchar *s1 = seq1ptr;

      // !!! careful: values written into pad cells may be wrong
      // !!! (overwrite each other, wrong scaling)

      if (!pssmPtr) {
	const uchar *s2 = seq2ptr;

	for (int i = 0; i < numCells; ++i) {
	  const double matchProb = substitutionProbs[*s1][*s2];
	  const double yM = bM0[i];
	  const double yD = bD0[i];
	  const double yI = bI0[i];
	  double ySum = yM * matchProb + yD * delInit + yI * insInit;

	  if (!globality || matchProb <= 0) ySum += scaledUnit;
	  // xxx matchProb should be 0 only at delimiters, but will be
	  // 0 for non-delimiters with severe mismatch scores

	  bM2[i] = ySum;
	  bD1[i] = ySum + yD * delNext;
	  bI1[i] = ySum + yI * insNext;

	  *mDout += fD0[i] * yD;
	  *mIout += fI0[i] * yI;

	  mDout++; mIout--;
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}
      } else {
	const ExpMatrixRow *p2 = pssmPtr;

	for (int i = 0; i < numCells; ++i) {
	  const double matchProb = (*p2)[*s1];
	  const double yM = bM0[i];
	  const double yD = bD0[i];
	  const double yI = bI0[i];
	  double ySum = yM * matchProb + yD * delInit + yI * insInit;

	  if (!globality || matchProb <= 0) ySum += scaledUnit;  // xxx

	  bM2[i] = ySum;
	  bD1[i] = ySum + yD * delNext;
	  bI1[i] = ySum + yI * insNext;

	  *mDout += fD0[i] * yD;
	  *mIout += fI0[i] * yI;

	  mDout++; mIout--;
	  s1 += seqIncrement;
	  p2 -= seqIncrement;
	}
      }

      if (antidiagonal == 0) break;

      oldPos = newPos;

      if ((antidiagonal + 2) % rescaleStep == 0 &&
	  antidiagonal + 2 < numAntidiagonals) {
	const double scale = rescales[antidiagonal / rescaleStep];
	rescaleBckProbs(diagPos, newPos, scale);
	scaledUnit *= scale;
      }

      --antidiagonal;
      const size_t newSeq1beg = xa.seq1start(antidiagonal);
      if (newSeq1beg < seq1beg) {
	seq1beg = newSeq1beg;
	seq1ptr -= seqIncrement;
      } else {
	seq2ptr -= seqIncrement;
	if (pssmPtr) pssmPtr -= seqIncrement;
      }
    }
  }

  double Centroid::dp_centroid( double gamma ){
    for( size_t k = 1; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t scoreEnd = xa.scoreEndIndex( k );
      double* X0 = &X[ scoreEnd ];
      size_t seq1pos = xa.seq1start( k );

      const double* const x0end = X0 + xa.numCellsAndPads( k );
      const size_t h = xa.hori( k, seq1pos );
      const size_t d = xa.diag( k, seq1pos );
      const double* X1 = &X[h];
      const double* X2 = &X[d];
      const double* fM2 = &fM[d];
      const double* bM2 = &bM[d];

      for (int i = 0; i < xdropPadLen; ++i) {
	*X0++ = -DINF;
      }

      do{
	const double matchProb = (*fM2++) * (*bM2++);
	const double s = ( gamma + 1 ) * matchProb - 1;
	const double oldX1 = *X1++;  // Added by MCF
	const double score = std::max( std::max( oldX1, *X1 ), *X2++ + s );
	//assert ( score >= 0 );
	updateScore ( score, k, seq1pos );
	*X0++ = score;
	seq1pos++;
      }while( X0 != x0end );
    }
    return bestScore;
  }

  void Centroid::traceback_centroid( std::vector< SegmentPair >& chunks,
				     double gamma ) const{
    //std::cout << "[c] bestAntiDiagonal=" << bestAntiDiagonal << ": bestPos1=" << bestPos1 << std::endl;

    size_t k = bestAntiDiagonal;
    size_t i = bestPos1;
    size_t oldPos1 = i;

    while( k > 0 ){
      const size_t h = xa.hori( k, i );
      const size_t v = xa.vert( k, i );
      const size_t d = xa.diag( k, i );
      const double matchProb = fM[d] * bM[d];
      const int m = maxIndex( X[d] + (gamma + 1) * matchProb - 1, X[h], X[v] );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }
      if( (m > 0 && oldPos1 != i) || k == 0 ){
	chunks.push_back( SegmentPair( i, k - i, oldPos1 - i ) );
      }
      if( m > 0 ){
	k -= 1;
	i -= (m == 1);
	oldPos1 = i;
      }
    }
  }

  double Centroid::dp_ama( double gamma ){
    mX1.assign ( numAntidiagonals + 2, 1.0 );
    mX2.assign ( numAntidiagonals + 2, 1.0 );

    for (size_t k = 0; k < numAntidiagonals; ++k) {
      size_t seq1pos = xa.seq1start(k);
      size_t seq2pos = k - seq1pos;
      size_t loopBeg = xa.diag(k, seq1pos);
      size_t loopEnd = loopBeg + xa.numCellsAndPads(k) - xdropPadLen;
      for (size_t i = loopBeg; i < loopEnd; ++i) {
	const double matchProb = fM[i] * bM[i];
	mX1[seq1pos++] -= matchProb;
	mX2[seq2pos--] -= matchProb;
      }
    }

    for( size_t k = 1; k < numAntidiagonals; ++k ){  // loop over antidiagonals
      const size_t scoreEnd = xa.scoreEndIndex( k );
      double* X0 = &X[ scoreEnd ];
      size_t seq1pos = xa.seq1start( k );
      size_t seq2pos = k - seq1pos;

      const double* const x0end = X0 + xa.numCellsAndPads( k );
      const size_t h = xa.hori( k, seq1pos );
      const size_t d = xa.diag( k, seq1pos );
      const double* X1 = &X[h];
      const double* X2 = &X[d];
      const double* fM2 = &fM[d];
      const double* bM2 = &bM[d];

      for (int i = 0; i < xdropPadLen; ++i) {
	*X0++ = -DINF;
      }

      do{
	const double matchProb = (*fM2++) * (*bM2++);
	const double thisD = mD[seq1pos];
	const double thisI = mI[seq2pos];
	const double thisXD = mX1[seq1pos] - thisD;
	const double thisXI = mX2[seq2pos] - thisI;
	const double s = 2 * gamma * matchProb - (thisXD + thisXI);
	const double u = gamma * thisD - thisXD;
	const double t = gamma * thisI - thisXI;
	const double oldX1 = *X1++;  // Added by MCF
	const double score = std::max(std::max(oldX1 + u, *X1 + t), *X2++ + s);
	updateScore ( score, k, seq1pos );
	*X0++ = score;
	seq1pos++;
	seq2pos--;
      }while( X0 != x0end );
    }

    return bestScore;
  }

  void Centroid::traceback_ama( std::vector< SegmentPair >& chunks,
			    double gamma ) const{
    //std::cout << "[c] bestAntiDiagonal=" << bestAntiDiagonal << ": bestPos1=" << bestPos1 << std::endl;

    size_t k = bestAntiDiagonal;
    size_t i = bestPos1;
    size_t oldPos1 = i;

    while( k > 0 ){
      const size_t j = k - i;
      const size_t h = xa.hori( k, i );
      const size_t v = xa.vert( k, i );
      const size_t d = xa.diag( k, i );
      const double matchProb = fM[d] * bM[d];
      const double thisD = mD[i];
      const double thisI = mI[j];
      const double thisXD = mX1[i] - thisD;
      const double thisXI = mX2[j] - thisI;
      const double s = 2 * gamma * matchProb - (thisXD + thisXI);
      const double t = gamma * thisI - thisXI;
      const double u = gamma * thisD - thisXD;
      const int m = maxIndex( X[d] + s, X[h] + u, X[v] + t );
      if( m == 0 ){
	k -= 2;
	i -= 1;
      }
      if( (m > 0 && oldPos1 != i) || k == 0 ){
	chunks.push_back( SegmentPair( i, k - i, oldPos1 - i ) );
      }
      if( m > 0 ){
	k -= 1;
	i -= (m == 1);
	oldPos1 = i;
      }
    }
  }

  void Centroid::getMatchAmbiguities(std::vector<char>& ambiguityCodes,
				     size_t seq1end, size_t seq2end,
				     size_t length) const {
    while (length) {
      size_t d = xa.diag(seq1end + seq2end, seq1end);
      double p = fM[d] * bM[d];
      ambiguityCodes.push_back(asciiProbability(p));
      --seq1end;  --seq2end;  --length;
    }
  }

  void Centroid::getDeleteAmbiguities(std::vector<char>& ambiguityCodes,
				      size_t seq1end, size_t seq1beg) const {
    for (size_t i = seq1end; i > seq1beg; --i)
      ambiguityCodes.push_back(asciiProbability(mD[i]));
  }

  void Centroid::getInsertAmbiguities(std::vector<char>& ambiguityCodes,
				      size_t seq2end, size_t seq2beg) const {
    for (size_t i = seq2end; i > seq2beg; --i)
      ambiguityCodes.push_back(asciiProbability(mI[i]));
  }

  static void countUncertainLetters(double *counts, double alignProb,
				    unsigned alphabetSize,
				    const double *probRatios,
				    const double *letterProbs) {
    double ratioParts[scoreMatrixRowSize];
    double sum = 0;
    for (unsigned letter = 0; letter < alphabetSize; ++letter) {
      double r = probRatios[letter] * letterProbs[letter];
      ratioParts[letter] = r;
      sum += r;
    }
    if (sum > 0) {
      const double mul = alignProb / sum;
      for (unsigned letter = 0; letter < alphabetSize; ++letter) {
	counts[letter] += mul * ratioParts[letter];
      }
    }
  }

  void Centroid::addExpectedCounts(size_t start2, bool isExtendFwd,
				   const const_dbl_ptr *substitutionProbs,
				   const GapCosts &gapCosts,
				   unsigned alphabetSize,
				   const dbl_ptr *substitutionCounts,
				   double *transitionCounts) {
    const double *letterProbs = 0;
    if (!letterProbsPerPosition.empty()) {
      letterProbs = &letterProbsPerPosition[0] + start2 * alphabetSize;
    }

    const int seqIncrement = isExtendFwd ? 1 : -1;
    int alphabetSizeIncrement = alphabetSize;
    if (!isExtendFwd) alphabetSizeIncrement *= -1;

    size_t antidiagonal = 0;
    size_t seq1beg = 0;
    size_t vertPos = xdropPadLen * 2;
    size_t thisPos = xdropPadLen * 3;

    double alignedLetterPairCount = 0;
    double delNextCount = 0;
    double insNextCount = 0;

    while (1) {
      const double *bM0 = &bM[thisPos];
      const double *bD0 = &bD[thisPos];
      const double *bI0 = &bI[thisPos];
      const double *fM0 = &fM[thisPos];
      const double *fD1 = &fD[vertPos - 1];
      const double *fI1 = &fI[vertPos];

      double dNextCount = 0;
      double iNextCount = 0;

      ++antidiagonal;
      const size_t nextPos = xa.scoreEndIndex(antidiagonal);
      const int numCells = nextPos - thisPos;
      const uchar *s1 = seq1ptr;

      if (!letterProbs) {
	const uchar *s2 = seq2ptr;
	for (int i = 0; i < numCells; ++i) {
	  const double alignProb = fM0[i] * bM0[i];
	  substitutionCounts[*s1][*s2] += alignProb;
	  alignedLetterPairCount += alignProb;
	  dNextCount += fD1[i] * bD0[i];
	  iNextCount += fI1[i] * bI0[i];
	  s1 += seqIncrement;
	  s2 -= seqIncrement;
	}
      } else {
	const double *lp2 = letterProbs;
	for (int i = 0; i < numCells; ++i) {
	  const double alignProb = fM0[i] * bM0[i];
	  const unsigned letter1 = *s1;
	  countUncertainLetters(substitutionCounts[letter1], alignProb,
				alphabetSize, substitutionProbs[letter1], lp2);
	  alignedLetterPairCount += alignProb;
	  dNextCount += fD1[i] * bD0[i];
	  iNextCount += fI1[i] * bI0[i];
	  s1 += seqIncrement;
	  lp2 -= alphabetSizeIncrement;
	}
      }

      if ((antidiagonal + 1) % rescaleStep == 0 &&
	  antidiagonal + 1 < numAntidiagonals) {
	const double mul = rescales[antidiagonal / rescaleStep];
	dNextCount *= mul;
	iNextCount *= mul;
      }

      delNextCount += dNextCount;
      insNextCount += iNextCount;

      if (antidiagonal == numAntidiagonals) break;

      vertPos = thisPos;
      thisPos = nextPos + xdropPadLen;

      const size_t newSeq1beg = xa.seq1start(antidiagonal);
      if (newSeq1beg > seq1beg) {
	seq1beg = newSeq1beg;
	seq1ptr += seqIncrement;
	++vertPos;
      } else {
	seq2ptr += seqIncrement;
	if (letterProbs) letterProbs += alphabetSizeIncrement;
      }
    }

    double delCount = 0;
    double insCount = 0;
    for (size_t i = 0; i < numAntidiagonals + 2; ++i) {
      delCount += mD[i];
      insCount += mI[i];
    }

    delNextCount *= gapCosts.delProbPieces[0].growProb;
    insNextCount *= gapCosts.insProbPieces[0].growProb;

    transitionCounts[0] += alignedLetterPairCount;
    transitionCounts[1] += delCount;  // deleted letter count
    transitionCounts[2] += insCount;  // inserted letter count
    transitionCounts[3] += delCount - delNextCount;  // delete open/close count
    transitionCounts[4] += insCount - insNextCount;  // insert open/close count
  }

}  // end namespace cbrc
