// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

#include "Alignment.hh"
#include "GeneticCode.hh"
#include "LastEvaluer.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"

#include <assert.h>

#include <algorithm>
#include <cstdio>  // sprintf

using namespace cbrc;

const char aminoTriplets[] =
  "Ala" "ala"
  "Cys" "cys"
  "Asp" "asp"
  "Glu" "glu"
  "Phe" "phe"
  "Gly" "gly"
  "His" "his"
  "Ile" "ile"
  "Lys" "lys"
  "Leu" "leu"
  "Met" "met"
  "Asn" "asn"
  "Pro" "pro"
  "Gln" "gln"
  "Arg" "arg"
  "Ser" "ser"
  "Thr" "thr"
  "Val" "val"
  "Trp" "trp"
  "Tyr" "tyr"
  "***" "***"
  "Xaa" "xaa";

// This writes a "size_t" integer into a char buffer ending at "end".
// It writes backwards from the end, because that's easier & faster.
static char *writeSize(char *end, size_t x) {
  do {
    --end;
    *end = '0' + x % 10;
    x /= 10;
  } while (x);
  return end;
}

// write x - y as a signed integer
static char *writeSignedDifference(char *end, size_t x, size_t y) {
  if (x >= y) {
    end = writeSize(end, x - y);
  } else {
    end = writeSize(end, y - x);
    *--end = '-';
  }
  return end;
}

class IntText {  // a text representation of an integer
public:
  IntText() {}
  explicit IntText(size_t x) { set(x); }
  void set(size_t x) { char *e = b + sizeof b; s = e - writeSize(e, x); }
  const char *begin() const { return b + sizeof b - s; }
  size_t size() const { return s; }
private:
  char b[31];
  unsigned char s;
};

class FloatText {  // a text representation of a floating-point number
public:
  FloatText() {}
  FloatText(const char *format, double x) { set(format, x); }
  void set(const char *format, double x) { s = std::sprintf(b, format, x); }
  const char *begin() const { return b; }
  size_t size() const { return s; }
private:
  char b[31];
  unsigned char s;
};

class Writer {  // writes characters to an output buffer
public:
  explicit Writer(char *startOfOutput) : p(startOfOutput) {}
  char *pointer() const { return p; }
  void fill(size_t n, char c) {  // write n copies of character c
    std::memset(p, c, n);
    p += n;
  }
  void copy(const char *from, size_t size) {
    std::memcpy(p, from, size);
    p += size;
  }
  Writer &operator<<(char c) {
    *p++ = c;
    return *this;
  }
  Writer &operator<<(const IntText &t) {
    copy(t.begin(), t.size());
    return *this;
  }
  Writer &operator<<(const FloatText &t) {
    copy(t.begin(), t.size());
    return *this;
  }
  Writer &operator<<(const std::string &t) {
    copy(t.c_str(), t.size());
    return *this;
  }
private:
  char *p;
};

AlignmentText Alignment::write(const MultiSequence& seq1,
			       const MultiSequence& seq2,
			       size_t seqNum2, const uchar* seqData2,
			       const Alphabet& alph, const Alphabet& dnaAlph,
			       int translationType, const uchar *codonToAmino,
			       const LastEvaluer& evaluer, int format,
			       const AlignmentExtras& extras) const {
  assert(!blocks.empty());

  if (format == 'm')
    return writeMaf(seq1, seq2, seqNum2, seqData2,
		    alph, dnaAlph, translationType, evaluer, extras);
  if (format == 't')
    return writeTab(seq1, seq2, seqNum2, translationType, evaluer, extras);
  else
    return writeBlastTab(seq1, seq2, seqNum2, seqData2, alph, translationType,
			 codonToAmino, evaluer, extras, format == 'B');
}

static size_t alignedColumnCount(const std::vector<SegmentPair> &blocks) {
  size_t c = 0;
  for (size_t i = 0; i < blocks.size(); ++i)
    c += blocks[i].size;
  return c;
}

static size_t matchCount(const std::vector<SegmentPair> &blocks,
			 const uchar *seq1, const uchar *seq2,
			 const uchar *map1, const uchar *map2) {
  // no special treatment of ambiguous bases/residues: same as NCBI BLAST
  size_t matches = 0;
  for (size_t i = 0; i < blocks.size(); ++i) {
    const SegmentPair &b = blocks[i];
    const uchar *x = seq1 + b.beg1();
    const uchar *y = seq2 + b.beg2();
    for (size_t j = 0; j < b.size; ++j)
      if (map1[x[j]] == map2[y[j]])
	++matches;
  }
  return matches;
}

static size_t writeBlocks(std::vector<char> &text,
			  const std::vector<SegmentPair> &blocks,
			  size_t frameSize2) {
  size_t s = blocks.size();
  text.resize(32 * 3 * s);
  char *end = &text[0] + text.size();
  char *e = end;
  for (size_t i = s; i --> 0; ) {
    const SegmentPair &y = blocks[i];
    if (y.size) e = writeSize(e, y.size);
    if (i > 0) {  // between each pair of aligned blocks:
      const SegmentPair &x = blocks[i - 1];
      if (y.size) *--e = ',';
      size_t gapBeg2 = aaToDna(x.end2(), frameSize2);
      size_t gapEnd2 = aaToDna(y.beg2(), frameSize2);
      e = writeSignedDifference(e, gapEnd2, gapBeg2);  // allow -1 frameshift
      *--e = ':';
      size_t gapBeg1 = x.end1();
      size_t gapEnd1 = y.beg1();
      e = writeSignedDifference(e, gapEnd1, gapBeg1);  // allow -1 frameshift
      if (x.size) *--e = ',';
    }
  }
  return end - e;
}

static const char *scoreFormat(bool isIntegerScores) {
  return isIntegerScores ? "%.0f" : "%.1f";
}

static char* writeTags( const LastEvaluer& evaluer, double queryLength,
			double score, double fullScore,
			char separator, char* out ){
  if( evaluer.isGood() ){
    double epa = evaluer.evaluePerArea( score );
    double area = evaluer.area( score, queryLength );
    *out++ = separator;
    out += std::sprintf( out, "EG2=%.2g", 1e18 * epa );
    *out++ = separator;
    out += std::sprintf( out, "E=%.2g", area * epa );
  }
  if( fullScore > 0 ){
    *out++ = separator;
    out += std::sprintf( out, "fullScore=%.3g", fullScore );
  }
  *out++ = '\n';
  return out;
}

AlignmentText Alignment::writeTab(const MultiSequence& seq1,
				  const MultiSequence& seq2,
				  size_t seqNum2, bool isTranslated,
				  const LastEvaluer& evaluer,
				  const AlignmentExtras& extras) const {
  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t seqNum1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(seqNum1);
  size_t seqLen1 = seq1.seqLen(seqNum1);

  size_t size2 = seq2.padLen(seqNum2);
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t seqStart2 = seq2.seqBeg(seqNum2) - seq2.padBeg(seqNum2);
  size_t seqLen2 = seq2.seqLen(seqNum2);

  FloatText sc(scoreFormat(extras.fullScore >= 0), score);
  std::string n1 = seq1.seqName(seqNum1);
  char strand1 = seq1.strand(seqNum1);
  std::string n2 = seq2.seqName(seqNum2);
  char strand2 = seq2.strand(seqNum2);
  IntText b1(alnBeg1 - seqStart1);
  IntText b2(alnBeg2 - seqStart2);
  IntText r1(alnEnd1 - alnBeg1);
  IntText r2(alnEnd2 - alnBeg2);
  IntText s1(seqLen1);
  IntText s2(seqLen2);

  std::vector<char> blockText;
  size_t blockLen = writeBlocks(blockText, blocks, frameSize2);

  const char t = '\t';
  char tags[256];
  char *tagEnd = writeTags(evaluer, seqLen2, score, extras.fullScore, t, tags);
  size_t tagLen = tagEnd - tags;

  size_t textLen = sc.size() + 1 +
    n1.size() + b1.size() + r1.size() + 1 + s1.size() + 5 +
    n2.size() + b2.size() + r2.size() + 1 + s2.size() + 5 + blockLen + tagLen;

  char *text = new char[textLen + 1];
  Writer w(text);
  w << sc << t;
  w << n1 << t << b1 << t << r1 << t << strand1 << t << s1 << t;
  w << n2 << t << b2 << t << r2 << t << strand2 << t << s2 << t;
  w.copy(&blockText[0] + blockText.size() - blockLen, blockLen);
  w.copy(tags, tagLen);
  w << '\0';

  return AlignmentText(seqNum2, alnBeg2, alnEnd2, strand2, score, 0, 0, text);
}

static void putLeft(Writer &w, const std::string &t, size_t width) {
  w << t;
  w.fill(width - t.size(), ' ');
  w << ' ';
}

static void putRight(Writer &w, const IntText &t, size_t width) {
  w.fill(width - t.size(), ' ');
  w << t;
  w << ' ';
}

// Write an "a" line
static char *writeMafLineA(char *out, double score, const LastEvaluer& evaluer,
			   double queryLength, double fullScore) {
  out += std::sprintf(out, "a score=");
  out += std::sprintf(out, scoreFormat(fullScore >= 0), score);
  return writeTags(evaluer, queryLength, score, fullScore, ' ', out);
}

// Write the first part of an "s" line:
static char *writeMafHeadS(char *out,
			   const std::string &n, size_t nw,
			   const IntText &b, size_t bw,
			   const IntText &r, size_t rw,
			   char strand,
			   const IntText &s, size_t sw) {
  Writer w(out);
  w << 's' << ' ';
  putLeft(w, n, nw);
  putRight(w, b, bw);
  putRight(w, r, rw);
  w << strand << ' ';
  putRight(w, s, sw);
  return w.pointer();
}

// Write the first part of a "q" line:
static char *writeMafHeadQ(char *out,
			   const std::string &n, size_t nw,
			   size_t qLineBlankLen) {
  Writer w(out);
  w << 'q' << ' ';
  putLeft(w, n, nw);
  w.fill(qLineBlankLen, ' ');
  return w.pointer();
}

// Write a "c" line
static void writeMafLineC(std::vector<char> &cLine,
			  const std::vector<double> &counts,
			  const Alphabet &alph, bool isCodon) {
  if (counts.empty()) return;
  unsigned alphSize2 = isCodon ? 64 : alph.size;
  unsigned numOfSubstitutions = scoreMatrixRowSize * scoreMatrixRowSize;
  size_t numOfTransitions = counts.size() - numOfSubstitutions;
  size_t numOfCounts = alph.size * alphSize2 + numOfTransitions;
  cLine.resize(2 + 32 * numOfCounts);
  char *e = &cLine[0];
  *e++ = 'c';

  for (unsigned i = 0; i < alph.size; ++i) {
    unsigned x = alph.numbersToLowercase[i];
    for (unsigned j = 0; j < alphSize2; ++j) {
      unsigned y = isCodon ? j : alph.numbersToLowercase[j];
      double c = counts[i * scoreMatrixRowSize + j];
      if (x != i) c += counts[x * scoreMatrixRowSize + j];
      if (y != j) c += counts[i * scoreMatrixRowSize + y];
      if (x != i && y != j) c += counts[x * scoreMatrixRowSize + y];
      e += std::sprintf(e, " %.3g", c);
    }
  }

  for (size_t i = 0; i < numOfTransitions; ++i) {
    e += std::sprintf(e, " %.3g", counts[numOfSubstitutions + i]);
  }

  *e++ = '\n';
  cLine.resize(e - &cLine[0]);
}

AlignmentText Alignment::writeMaf(const MultiSequence& seq1,
				  const MultiSequence& seq2,
				  size_t seqNum2, const uchar* seqData2,
				  const Alphabet& alph,
				  const Alphabet& dnaAlph,
				  int translationType,
				  const LastEvaluer& evaluer,
				  const AlignmentExtras& extras) const {
  bool isCodon = (translationType == 2);
  double fullScore = extras.fullScore;
  const std::vector<char>& columnProbSymbols = extras.columnAmbiguityCodes;

  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t seqNum1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(seqNum1);
  size_t seqLen1 = seq1.seqLen(seqNum1);

  size_t size2 = seq2.padLen(seqNum2);
  size_t frameSize2 = translationType ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t seqOrigin2 = seq2.padBeg(seqNum2);
  size_t seqStart2 = seq2.seqBeg(seqNum2) - seqOrigin2;
  size_t seqLen2 = seq2.seqLen(seqNum2);

  char aLine[256];
  char *aLineEnd = writeMafLineA(aLine, score, evaluer, seqLen2, fullScore);
  size_t aLineLen = aLineEnd - aLine;

  const std::string n1 = seq1.seqName(seqNum1);
  char strand1 = seq1.strand(seqNum1);
  const std::string n2 = seq2.seqName(seqNum2);
  char strand2 = seq2.strand(seqNum2);
  IntText b1(alnBeg1 - seqStart1);
  IntText b2(alnBeg2 - seqStart2);
  IntText r1(alnEnd1 - alnBeg1);
  IntText r2(alnEnd2 - alnBeg2);
  IntText s1(seqLen1);
  IntText s2(seqLen2);

  const size_t nw = std::max( n1.size(), n2.size() );
  const size_t bw = std::max( b1.size(), b2.size() );
  const size_t rw = std::max( r1.size(), r2.size() );
  const size_t sw = std::max( s1.size(), s2.size() );

  size_t qLineBlankLen = bw + 1 + rw + 3 + sw + 1;
  size_t pLineBlankLen = nw + 1 + qLineBlankLen;
  size_t sLineLen = 2 + pLineBlankLen + numColumns(frameSize2, isCodon) + 1;

  std::vector<char> cLine;
  writeMafLineC(cLine, extras.expectedCounts, alph, isCodon);

  size_t qualsPerBase1 = seq1.qualsPerLetter();
  size_t qualsPerBase2 = seq2.qualsPerLetter();
  bool isQuals1 = qualsPerBase1 && (translationType != 2);
  bool isQuals2 = qualsPerBase2 && (translationType != 1);

  size_t sLineNum = 2 + isQuals1 + isQuals2 + !columnProbSymbols.empty();
  size_t textLen = aLineLen + sLineLen * sLineNum + cLine.size() + 1;
  char *text = new char[textLen + 1];

  char *dest = std::copy(aLine, aLineEnd, text);

  dest = writeMafHeadS(dest, n1, nw, b1, bw, r1, rw, strand1, s1, sw);
  dest = writeTopSeq(dest, seq1.seqReader(), alph, 0, frameSize2, isCodon);
  *dest++ = '\n';

  if (isQuals1) {
    dest = writeMafHeadQ(dest, n1, nw, qLineBlankLen);
    const uchar *q = seq1.qualityReader();
    dest = writeTopSeq(dest, q, alph, qualsPerBase1, frameSize2, isCodon);
    *dest++ = '\n';
  }

  dest = writeMafHeadS(dest, n2, nw, b2, bw, r2, rw, strand2, s2, sw);
  if (isCodon) seqData2 = seq2.seqReader() + seqOrigin2;
  const Alphabet &alph2 = isCodon ? dnaAlph : alph;
  dest = writeBotSeq(dest, seqData2, alph2, 0, frameSize2, isCodon);
  *dest++ = '\n';

  if (isQuals2) {
    dest = writeMafHeadQ(dest, n2, nw, qLineBlankLen);
    const uchar *q = seq2.qualityReader() + seqOrigin2 * qualsPerBase2;
    dest = writeBotSeq(dest, q, alph2, qualsPerBase2, frameSize2, isCodon);
    *dest++ = '\n';
  }

  if (!columnProbSymbols.empty()) {
    Writer w(dest);
    w << 'p' << ' ';
    w.fill(pLineBlankLen, ' ');
    dest = w.pointer();
    dest = writeColumnProbs(dest, &columnProbSymbols[0], frameSize2, isCodon);
    *dest++ = '\n';
  }

  Writer w(dest);
  if (!cLine.empty()) w.copy(&cLine[0], cLine.size());
  w << '\n' << '\0';  // blank line afterwards

  return AlignmentText(seqNum2, alnBeg2, alnEnd2, strand2, score, 0, 0, text);
}

AlignmentText Alignment::writeBlastTab(const MultiSequence& seq1,
				       const MultiSequence& seq2,
				       size_t seqNum2, const uchar* seqData2,
				       const Alphabet& alph,
				       int translationType,
				       const uchar *codonToAmino,
				       const LastEvaluer& evaluer,
				       const AlignmentExtras& extras,
				       bool isExtraColumns) const {
  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t seqNum1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(seqNum1);
  size_t seqLen1 = seq1.seqLen(seqNum1);
  char strand1 = seq1.strand(seqNum1);

  size_t size2 = seq2.padLen(seqNum2);
  size_t frameSize2 = translationType ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t seqStart2 = seq2.seqBeg(seqNum2) - seq2.padBeg(seqNum2);
  size_t seqLen2 = seq2.seqLen(seqNum2);
  char strand2 = seq2.strand(seqNum2);

  size_t alnSize = numColumns(frameSize2, false);
  const uchar *map1 = alph.numbersToUppercase;
  const uchar *map2 = (translationType == 2) ? codonToAmino : map1;
  size_t matches = matchCount(blocks, seq1.seqReader(), seqData2, map1, map2);
  size_t mismatches = alignedColumnCount(blocks) - matches;
  size_t gapOpens = blocks.size() - 1;  // xxx ???
  double matchPercent = 100.0 * matches / alnSize;

  size_t blastAlnBeg1 = alnBeg1 + 1;  // 1-based coordinate
  size_t blastAlnEnd1 = alnEnd1;
  if (strand1 == '-') {
    blastAlnBeg1 = seqStart1 + seqLen1 - alnBeg1;
    blastAlnEnd1 = seqStart1 + seqLen1 - alnEnd1 + 1;  // 1-based coordinate
    seqStart1 = 0;
  }
  size_t blastAlnBeg2 = alnBeg2 + 1;  // 1-based coordinate
  size_t blastAlnEnd2 = alnEnd2;
  if (strand2 == '-') {
    blastAlnBeg2 = size2 - alnBeg2;
    blastAlnEnd2 = size2 - alnEnd2 + 1;  // 1-based coordinate
    /*
    if (!isTranslated) {  // xxx this makes it more like BLAST
      std::swap(blastAlnBeg1, blastAlnEnd1);
      std::swap(blastAlnBeg2, blastAlnEnd2);
    }
    */
  }

  std::string n1 = seq1.seqName(seqNum1);
  std::string n2 = seq2.seqName(seqNum2);
  FloatText mp("%.2f", matchPercent);
  IntText as(alnSize);
  IntText mm(mismatches);
  IntText go(gapOpens);
  IntText b1(blastAlnBeg1 - seqStart1);
  IntText b2(blastAlnBeg2 - seqStart2);
  IntText e1(blastAlnEnd1 - seqStart1);
  IntText e2(blastAlnEnd2 - seqStart2);
  FloatText ev;
  FloatText bs;
  if( evaluer.isGood() ){
    double area = evaluer.area( score, seqLen2 );
    double epa = evaluer.evaluePerArea( score );
    double bitScore = evaluer.bitScore( score );
    ev.set("%.2g", area * epa);
    bs.set("%.3g", bitScore);
  }
  IntText s1(seqLen1);
  IntText s2(seqLen2);
  FloatText sc;

  size_t s =
    n2.size() + n1.size() + mp.size() + as.size() + mm.size() + go.size() +
    b2.size() + e2.size() + b1.size() + e1.size() + 10;
  if (evaluer.isGood()) s += ev.size() + bs.size() + 2;
  if (isExtraColumns) {
    sc.set(scoreFormat(extras.fullScore >= 0), score);
    s += s1.size() + s2.size() + sc.size() + 3;
  }

  char *text = new char[s + 1];
  Writer w(text);
  const char t = '\t';
  w << n2 << t << n1 << t << mp << t << as << t << mm << t << go << t
    << b2 << t << e2 << t << b1 << t << e1;
  if (evaluer.isGood()) w << t << ev << t << bs;
  if (isExtraColumns)   w << t << s2 << t << s1 << t << sc;
  w << '\n' << '\0';

  return AlignmentText(seqNum2, alnBeg2, alnEnd2, strand2, score,
		       alnSize, matches, text);
}

static void setFrameCoords(size_t &seqBeg, size_t &seqEnd, size_t &frameshift,
			   size_t alnBeg, size_t alnEnd,
			   size_t frameSize, bool isCodon) {
  seqBeg = alnBeg;
  seqEnd = alnEnd;
  frameshift = 0;
  if (frameSize) {
    seqBeg = aaToDna(alnBeg, frameSize);
    seqEnd = aaToDna(alnEnd, frameSize);
    if (seqBeg > seqEnd) frameshift = 2;  // weird reverse frameshift
    if (!isCodon) {
      size_t len = seqEnd - seqBeg;
      seqBeg = alnEnd - (len + 1) / 3;
      seqEnd = alnEnd;
      frameshift = (len + 3) % 3;
    }
  }
}

size_t Alignment::numColumns(size_t frameSize, bool isCodon) const {
  const int aaLen = isCodon ? 3 : 1;
  size_t num = 0;

  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // length of unaligned chunk of top sequence (gaps in bottom sequence):
      num += (y.beg1() - x.end1()) * aaLen;

      // length of unaligned chunk of bottom sequence (gaps in top sequence):
      size_t beg, end, shift;
      setFrameCoords(beg, end, shift, x.end2(), y.beg2(), frameSize, isCodon);
      if (shift) ++num;
      if (beg < end) num += end - beg;
    }

    num += y.size * aaLen;  // length of aligned chunk
  }

  return num;
}

static char* writeGaps( char* dest, size_t num ){
  char* end = dest + num;
  while( dest < end ){
    *dest++ = '-';
  }
  return dest;
}

static char *writeQuals(char *dest, const uchar *qualities,
			size_t beg, size_t end, size_t qualsPerBase) {
  for( size_t i = beg; i < end; ++i ){
    const uchar* q = qualities + i * qualsPerBase;
    *dest++ = *std::max_element( q, q + qualsPerBase );
  }
  return dest;
}

static char *writeAmino(char *dest, const uchar *beg, const uchar *end,
			const uchar *toUnmasked) {
  for (const uchar *i = beg; i < end; ++i) {
    unsigned x = toUnmasked[*i];
    bool isMasked = (x != *i);
    if (x > 21) x = 21;
    unsigned y = (x * 2 + isMasked) * 3;
    *dest++ = aminoTriplets[y];
    *dest++ = aminoTriplets[y + 1];
    *dest++ = aminoTriplets[y + 2];
  }
  return dest;
}

static char *writeSeq(char *dest, const uchar *seq, size_t beg, size_t end,
		      const Alphabet &alph, size_t qualsPerBase, bool isAA) {
  return qualsPerBase ? writeQuals(dest, seq, beg, end, qualsPerBase)
    : isAA ? writeAmino(dest, seq + beg, seq + end, alph.numbersToUppercase)
    :        alph.rtCopy(seq + beg, seq + end, dest);
}

char *Alignment::writeTopSeq(char *dest, const uchar *seq,
			     const Alphabet &alph, size_t qualsPerBase,
			     size_t frameSize, bool isCodon) const {
  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // append unaligned chunk of top sequence:
      dest =
	writeSeq(dest, seq, x.end1(), y.beg1(), alph, qualsPerBase, isCodon);

      // append gaps for unaligned chunk of bottom sequence:
      size_t beg, end, shift;
      setFrameCoords(beg, end, shift, x.end2(), y.beg2(), frameSize, isCodon);
      if (shift) *dest++ = '-';
      if (beg < end) dest = writeGaps(dest, end - beg);
    }

    // append aligned chunk of top sequence:
    dest =
      writeSeq(dest, seq, y.beg1(), y.end1(), alph, qualsPerBase, isCodon);
  }

  return dest;
}

char *Alignment::writeBotSeq(char *dest, const uchar *seq,
			     const Alphabet &alph, size_t qualsPerBase,
			     size_t frameSize, bool isCodon) const {
  const int aaLen = isCodon ? 3 : 1;

  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // append gaps for unaligned chunk of top sequence:
      dest = writeGaps(dest, (y.beg1() - x.end1()) * aaLen);

      //append unaligned chunk of bottom sequence:
      size_t beg, end, shift;
      setFrameCoords(beg, end, shift, x.end2(), y.beg2(), frameSize, isCodon);
      if (shift == 1) *dest++ = '\\';
      if (shift == 2) *dest++ = '/';
      dest = writeSeq(dest, seq, beg, end, alph, qualsPerBase, false);
    }

    // append aligned chunk of bottom sequence:
    size_t beg = isCodon ? aaToDna(y.beg2(), frameSize) : y.beg2();
    size_t end = isCodon ? aaToDna(y.end2(), frameSize) : y.end2();
    dest = writeSeq(dest, seq, beg, end, alph, qualsPerBase, false);
  }

  return dest;
}

static char *writeProbSymbols(char *dest, const char *probSymbols,
			      size_t size, int aaLen) {
  for (size_t i = 0; i < size; ++i)
    for (int j = 0; j < aaLen; ++j)
      *dest++ = probSymbols[i];
  return dest;
}

char *Alignment::writeColumnProbs(char *dest, const char *probSymbols,
				  size_t frameSize, bool isCodon) const {
  const int aaLen = isCodon ? 3 : 1;

  for (size_t i = 0; i < blocks.size(); ++i) {
    const SegmentPair &y = blocks[i];
    if (i > 0) {  // between each pair of aligned blocks:
      const SegmentPair &x = blocks[i - 1];

      size_t topLen = y.beg1() - x.end1();
      dest = writeProbSymbols(dest, probSymbols, topLen, aaLen);
      probSymbols += topLen;

      size_t beg, end, shift;
      setFrameCoords(beg, end, shift, x.end2(), y.beg2(), frameSize, isCodon);
      size_t rawLen = end - beg;
      size_t botLen = (rawLen + aaLen / 3) / aaLen;
      if (shift) *dest++ = '-';
      dest = writeProbSymbols(dest, probSymbols, botLen, aaLen);
      probSymbols += botLen;
      if (beg < end && rawLen % aaLen == 1) *dest++ = '-';
      if (beg < end && rawLen % aaLen == 2) --dest;
    }

    dest = writeProbSymbols(dest, probSymbols, y.size, aaLen);
    probSymbols += y.size;
  }

  return dest;
}
