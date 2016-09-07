// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

#include "Alignment.hh"
#include "GeneticCode.hh"
#include "LastEvaluer.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"
#include <algorithm>
#include <cassert>
#include <cstdio>  // sprintf

using namespace cbrc;

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
			       size_t seqNum2, char strand,
			       const uchar* seqData2,
			       bool isTranslated, const Alphabet& alph,
			       const LastEvaluer& evaluer, int format,
			       const AlignmentExtras& extras) const {
  assert( !blocks.empty() );

  if( format == 'm' )
    return writeMaf( seq1, seq2, seqNum2, strand, seqData2,
		     isTranslated, alph, evaluer, extras );
  if( format == 't' )
    return writeTab( seq1, seq2, seqNum2, strand,
		     isTranslated, evaluer, extras );
  else
    return writeBlastTab( seq1, seq2, seqNum2, strand, seqData2,
			  isTranslated, alph, evaluer, format == 'B' );
}

static size_t alignedColumnCount(const std::vector<SegmentPair> &blocks) {
  size_t c = 0;
  for (size_t i = 0; i < blocks.size(); ++i)
    c += blocks[i].size;
  return c;
}

static size_t matchCount(const std::vector<SegmentPair> &blocks,
			 const uchar *seq1, const uchar *seq2,
			 const uchar *numbersToUppercase) {
  // no special treatment of ambiguous bases/residues: same as NCBI BLAST
  size_t matches = 0;
  for (size_t i = 0; i < blocks.size(); ++i) {
    const SegmentPair &b = blocks[i];
    const uchar *x = seq1 + b.beg1();
    const uchar *y = seq2 + b.beg2();
    for (size_t j = 0; j < b.size; ++j)
      if (numbersToUppercase[x[j]] == numbersToUppercase[y[j]])
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

static char* writeTags( const LastEvaluer& evaluer, double queryLength,
			int score, double fullScore,
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
				  size_t seqNum2, char strand,
				  bool isTranslated,
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

  IntText sc(score);
  std::string n1 = seq1.seqName(seqNum1);
  std::string n2 = seq2.seqName(seqNum2);
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
  w << n1 << t << b1 << t << r1 << t << '+'    << t << s1 << t;
  w << n2 << t << b2 << t << r2 << t << strand << t << s2 << t;
  w.copy(&blockText[0] + blockText.size() - blockLen, blockLen);
  w.copy(tags, tagLen);
  w << '\0';

  return AlignmentText(seqNum2, alnBeg2, alnEnd2, strand, score, 0, 0, text);
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
static char *writeMafLineA(char *out, int score, const LastEvaluer& evaluer,
			   double queryLength, double fullScore) {
  out += std::sprintf(out, "a score=%d", score);
  return writeTags(evaluer, queryLength, score, fullScore, ' ', out);
}

// Write the first part of an "s" line:
static void writeMafHeadS(char *out,
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
}

// Write the first part of a "q" line:
static void writeMafHeadQ(char *out,
			  const std::string &n, size_t nw,
			  size_t qLineBlankLen) {
  Writer w(out);
  w << 'q' << ' ';
  putLeft(w, n, nw);
  w.fill(qLineBlankLen, ' ');
}

// Write the first part of a "p" line:
static void writeMafHeadP(char *out, size_t pLineBlankLen) {
  Writer w(out);
  w << 'p' << ' ';
  w.fill(pLineBlankLen, ' ');
}

// Write a "c" line
static void writeMafLineC(std::vector<char> &cLine,
			  const std::vector<double> &counts) {
  size_t s = counts.size();
  if (s == 0) return;
  cLine.resize(2 + 32 * s);
  char *e = &cLine[0];
  *e++ = 'c';
  for (size_t i = 0; i < s; ++i)
    e += std::sprintf(e, " %.3g", counts[i]);
  *e++ = '\n';
  cLine.resize(e - &cLine[0]);
}

AlignmentText Alignment::writeMaf(const MultiSequence& seq1,
				  const MultiSequence& seq2,
				  size_t seqNum2, char strand,
				  const uchar* seqData2,
				  bool isTranslated, const Alphabet& alph,
				  const LastEvaluer& evaluer,
				  const AlignmentExtras& extras) const {
  double fullScore = extras.fullScore;
  const std::vector<uchar>& columnAmbiguityCodes = extras.columnAmbiguityCodes;

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

  char aLine[256];
  char *aLineEnd = writeMafLineA(aLine, score, evaluer, seqLen2, fullScore);
  size_t aLineLen = aLineEnd - aLine;

  const std::string n1 = seq1.seqName(seqNum1);
  const std::string n2 = seq2.seqName(seqNum2);
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
  size_t headLen = 2 + pLineBlankLen;
  size_t sLineLen = headLen + numColumns( frameSize2 ) + 1;

  std::vector<char> cLine;
  writeMafLineC(cLine, extras.expectedCounts);

  size_t qualsPerBase1 = seq1.qualsPerLetter();
  size_t qualsPerBase2 = seq2.qualsPerLetter();
  bool isQuals1 = qualsPerBase1;
  bool isQuals2 = qualsPerBase2 && !isTranslated;
  // for translated alignment: don't write untranslated quality data

  size_t sLineNum = 2 + isQuals1 + isQuals2 + !columnAmbiguityCodes.empty();
  size_t textLen = aLineLen + sLineLen * sLineNum + cLine.size() + 1;
  char *text = new char[textLen + 1];

  char *dest = std::copy(aLine, aLineEnd, text);

  writeMafHeadS(dest, n1, nw, b1, bw, r1, rw, '+', s1, sw);
  dest = writeTopSeq(seq1.seqReader(), alph, 0, frameSize2, dest + headLen);
  *dest++ = '\n';

  if (isQuals1) {
    writeMafHeadQ(dest, n1, nw, qLineBlankLen);
    const uchar *q = seq1.qualityReader();
    dest = writeTopSeq(q, alph, qualsPerBase1, frameSize2, dest + headLen);
    *dest++ = '\n';
  }

  writeMafHeadS(dest, n2, nw, b2, bw, r2, rw, strand, s2, sw);
  dest = writeBotSeq(seqData2, alph, 0, frameSize2, dest + headLen);
  *dest++ = '\n';

  if (isQuals2) {
    writeMafHeadQ(dest, n2, nw, qLineBlankLen);
    const uchar *q =
      seq2.qualityReader() + seq2.padBeg(seqNum2) * qualsPerBase2;
    dest = writeBotSeq(q, alph, qualsPerBase2, frameSize2, dest + headLen);
    *dest++ = '\n';
  }

  if (!columnAmbiguityCodes.empty()) {
    writeMafHeadP(dest, pLineBlankLen);
    dest = copy(columnAmbiguityCodes.begin(),
		columnAmbiguityCodes.end(), dest + headLen);
    *dest++ = '\n';
  }

  dest = copy(cLine.begin(), cLine.end(), dest);

  *dest++ = '\n';  // blank line afterwards
  *dest++ = '\0';

  return AlignmentText(seqNum2, alnBeg2, alnEnd2, strand, score, 0, 0, text);
}

AlignmentText Alignment::writeBlastTab(const MultiSequence& seq1,
				       const MultiSequence& seq2,
				       size_t seqNum2, char strand,
				       const uchar* seqData2,
				       bool isTranslated, const Alphabet& alph,
				       const LastEvaluer& evaluer,
				       bool isExtraColumns) const {
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

  size_t alnSize = numColumns( frameSize2 );
  size_t matches = matchCount( blocks, seq1.seqReader(), seqData2,
			       alph.numbersToUppercase );
  size_t mismatches = alignedColumnCount(blocks) - matches;
  size_t gapOpens = blocks.size() - 1;
  double matchPercent = 100.0 * matches / alnSize;

  size_t blastAlnBeg1 = alnBeg1 + 1;  // 1-based coordinate
  size_t blastAlnEnd1 = alnEnd1;
  size_t blastAlnBeg2 = alnBeg2 + 1;  // 1-based coordinate
  size_t blastAlnEnd2 = alnEnd2;
  if (strand == '-') {
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

  size_t s =
    n2.size() + n1.size() + mp.size() + as.size() + mm.size() + go.size() +
    b2.size() + e2.size() + b1.size() + e1.size() + 10;
  if (evaluer.isGood()) s += ev.size() + bs.size() + 2;
  if (isExtraColumns)   s += s1.size() + s2.size() + 2;

  char *text = new char[s + 1];
  Writer w(text);
  const char t = '\t';
  w << n2 << t << n1 << t << mp << t << as << t << mm << t << go << t
    << b2 << t << e2 << t << b1 << t << e1;
  if (evaluer.isGood()) w << t << ev << t << bs;
  if (isExtraColumns)   w << t << s2 << t << s1;
  w << '\n' << '\0';

  return AlignmentText(seqNum2, alnBeg2, alnEnd2, strand, score,
		       alnSize, matches, text);
}

size_t Alignment::numColumns( size_t frameSize ) const{
  size_t num = 0;

  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // length of unaligned chunk of top sequence (gaps in bottom sequence):
      num += y.beg1() - x.end1();

      // length of unaligned chunk of bottom sequence (gaps in top sequence):
      size_t gap2, frameshift2;
      sizeAndFrameshift( x.end2(), y.beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) ++num;
      num += gap2;
    }

    num += y.size;  // length of aligned chunk
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

static char* writeQuals( const uchar* qualities, size_t beg, size_t end,
			 size_t qualsPerBase, char* dest ){
  for( size_t i = beg; i < end; ++i ){
    const uchar* q = qualities + i * qualsPerBase;
    *dest++ = *std::max_element( q, q + qualsPerBase );
  }
  return dest;
}

static char* writeSeq( const uchar* seq, size_t beg, size_t end,
		       const Alphabet& alph, size_t qualsPerBase, char* dest ){
  return qualsPerBase ? writeQuals( seq, beg, end, qualsPerBase, dest )
    :                   alph.rtCopy( seq + beg, seq + end, dest );
}

char* Alignment::writeTopSeq( const uchar* seq, const Alphabet& alph,
			      size_t qualsPerBase, size_t frameSize,
			      char* dest ) const{
  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // append unaligned chunk of top sequence:
      dest = writeSeq( seq, x.end1(), y.beg1(), alph, qualsPerBase, dest );

      // append gaps for unaligned chunk of bottom sequence:
      size_t gap2, frameshift2;
      sizeAndFrameshift( x.end2(), y.beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) *dest++ = '-';
      dest = writeGaps( dest, gap2 );
    }

    // append aligned chunk of top sequence:
    dest = writeSeq( seq, y.beg1(), y.end1(), alph, qualsPerBase, dest);
  }

  return dest;
}

char* Alignment::writeBotSeq( const uchar* seq, const Alphabet& alph,
			      size_t qualsPerBase, size_t frameSize,
			      char* dest ) const{
  for( size_t i = 0; i < blocks.size(); ++i ){
    const SegmentPair& y = blocks[i];
    if( i > 0 ){  // between each pair of aligned blocks:
      const SegmentPair& x = blocks[i - 1];

      // append gaps for unaligned chunk of top sequence:
      dest = writeGaps( dest, y.beg1() - x.end1() );

      //append unaligned chunk of bottom sequence:
      size_t gap2, frameshift2;
      sizeAndFrameshift( x.end2(), y.beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 == 1 ) *dest++ = '\\';
      if( frameshift2 == 2 ) *dest++ = '/';
      size_t chunkBeg2 = y.beg2() - gap2;
      dest = writeSeq( seq, chunkBeg2, y.beg2(), alph, qualsPerBase, dest );
    }

    // append aligned chunk of bottom sequence:
    dest = writeSeq( seq, y.beg2(), y.end2(), alph, qualsPerBase, dest );
  }

  return dest;
}
