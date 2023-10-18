// Copyright 2013, 2014 Martin C. Frith

#include "cbrc_unsplit_alignment.hh"

#include <ctype.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <string.h>

#include <algorithm>
//#include <iostream>
#include <stdexcept>
#include <string>

static void err(const std::string& s) {
  throw std::runtime_error(s);
}

static bool isGraph(char c) {
  return c > ' ';  // faster than std::isgraph
}

static bool isSpace(char c) {
  return c > 0 && c <= ' ';  // faster than std::isspace
}

static bool isDigit(char c) {
  return c >= '0' && c <= '9';
}

namespace cbrc {
    
static const char *readUint(const char *c, unsigned &x) {
  if (!c) return 0;

  // faster than std::strtoul
  while (isSpace(*c)) ++c;
  if (!isDigit(*c)) return 0;
  unsigned z = *c++ - '0';
  while (isDigit(*c)) {
    if (z > UINT_MAX / 10) return 0;
    unsigned digit = *c++ - '0';
    z = z * 10 + digit;
    if (z < digit) return 0;
  }

  x = z;
  return c;
}

static const char *readChar(const char *c, char &d) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  if (*c == 0) return 0;
  d = *c++;
  return c;
}

static const char *skipWord(const char *c) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  const char *e = c;
  while (isGraph(*e)) ++e;
  if (e == c) return 0;
  return e;
}

static const char *skipSpace(const char *c) {
  if (!c) return 0;
  while (isSpace(*c)) ++c;
  return c;
}

static char *rskipSpace(char *c) {
  while (isSpace(*(c-1))) --c;
  return c;
}

struct Complement {
  unsigned char c[UCHAR_MAX+1];
  Complement() {
    static const char x[] = "ACGTRYSWKMBDHVN";
    static const char y[] = "TGCAYRSWMKVHDBN";
    for (unsigned i = 0; i < UCHAR_MAX+1; ++i) c[i] = i;
    for (unsigned i = 0; x[i] && y[i]; ++i) {
      c[toupper(x[i])] = toupper(y[i]);
      c[tolower(x[i])] = tolower(y[i]);
    }
  }
  unsigned char operator()(unsigned char x) {
    return c[x];
  }
};
static Complement complement;

static void reverseComplement(char *beg, char *end) {
  std::reverse(beg, end);
  std::transform(beg, end, beg, complement);
}

static void parseSeqLine(char *line, char *lineEnd, const char *&seqName,
			 unsigned &start, unsigned &span, char &strand,
			 unsigned &seqLen, char *&aln, char *&alnEnd) {
  seqName = skipSpace(skipWord(line));
  const char *nameEnd = skipWord(seqName);
  const char *s = nameEnd;
  s = readUint(s, start);
  s = readUint(s, span);
  s = readChar(s, strand);
  s = readUint(s, seqLen);
  s = skipSpace(s);
  if (!s || s >= lineEnd) err(std::string("bad MAF line: ") + line);
  aln = line + (s - line);
  alnEnd = lineEnd;
  line[nameEnd - line] = 0;
  *lineEnd = 0;  // trim any trailing whitespace
}

static void parseQualLine(char *line, char *lineEnd, char *&qual) {
  const char *s = skipSpace(skipWord(skipWord(line)));
  if (!s || s >= lineEnd) err(std::string("bad MAF line: ") + line);
  qual = line + (s - line);
  *lineEnd = 0;  // trim any trailing whitespace
}

void UnsplitAlignment::init(bool isTopSeqQuery) {
  const unsigned rankOfQrySeq = 2 - isTopSeqQuery;
  const unsigned rankOfRefSeq = isTopSeqQuery + 1;

  char refStrand, qryStrand;
  unsigned refSpan, qrySpan;
  unsigned refSeqLen, qrySeqLen;
  char *refAln;
  char *refAlnEnd;
  char *qryAln;
  char *qryAlnEnd;
  char *qual = 0;
  unsigned sLineCount = 0;

  for (char **i = linesBeg; i < linesEnd; ++i) {
    char *line = i[0];
    char *lineEnd = rskipSpace(i[1] - 1);
    if (line[0] == 's') {
      ++sLineCount;
      if (sLineCount == rankOfRefSeq) {
	parseSeqLine(line, lineEnd, rname, rstart, refSpan, refStrand,
		     refSeqLen, refAln, refAlnEnd);
      } else {
	parseSeqLine(line, lineEnd, qname, qstart, qrySpan, qryStrand,
		     qrySeqLen, qryAln, qryAlnEnd);
      }
    } else if (line[0] == 'q') {
      if (sLineCount != rankOfQrySeq) {
	err("I can only handle quality data for the query sequence");
      }
      parseQualLine(line, lineEnd, qual);
      if (qryStrand == '-') std::reverse(qual, lineEnd);
    }
  }

  if (sLineCount != 2) err("I need 2 sequences per alignment");

  if (qryStrand == '-') {
    rstart = refSeqLen - rstart - refSpan;
    qstart = qrySeqLen - qstart - qrySpan;
    reverseComplement(refAln, refAlnEnd);
    reverseComplement(qryAln, qryAlnEnd);
  }

  qstrand = (qryStrand != refStrand) * 2 + (qryStrand == '-');

  rend = rstart + refSpan;
  qend = qstart + qrySpan;
  ralign = refAln;
  qalign = qryAln;
  qQual = qual;
}

static unsigned seqPosFromAlnPos(unsigned alnPos, const char *aln) {
  return alnPos - std::count(aln, aln + alnPos, '-');
}

static unsigned nthBasePrefix(const char* sequenceWithGapsBeg, unsigned n) {
  for (unsigned i = 0; /* noop */; ++i)
    if (sequenceWithGapsBeg[i] != '-') {
      if (n > 0) --n;
      else return i;
    }
}

static unsigned nthBaseSuffix(const char *sequenceWithGapsEnd, unsigned n) {
  for (unsigned i = 0; /* noop */; ++i)
    if (*(sequenceWithGapsEnd - 1 - i) != '-') {
      if (n > 0) --n;
      else return i;
    }
}

void mafSliceBeg(const char* rAln, const char* qAln,
		 unsigned qBeg, unsigned& qSliceBeg, unsigned& alnBeg) {
  if (qSliceBeg < qBeg) {
    qSliceBeg = qBeg;
    alnBeg = 0;
  } else {
    alnBeg = nthBasePrefix(qAln, qSliceBeg - qBeg);
  }
  unsigned numInserts = nthBasePrefix(rAln + alnBeg, 0);
  alnBeg += numInserts;
  qSliceBeg += numInserts;
}

void mafSliceEnd(const char* rAln, const char* qAln,
		 unsigned qEnd, unsigned& qSliceEnd, unsigned& alnEnd) {
  unsigned alnLength = strlen(qAln);
  if (qSliceEnd > qEnd) {
    qSliceEnd = qEnd;
    alnEnd = alnLength;
  } else {
    alnEnd = alnLength - nthBaseSuffix(qAln + alnLength, qEnd - qSliceEnd);
  }
  unsigned numInserts = nthBaseSuffix(rAln + alnEnd, 0);
  alnEnd -= numInserts;
  qSliceEnd -= numInserts;
}

static const char *updateFieldWidth(const char *p, int &width) {
  while (isSpace(*p)) ++p;
  const char *b = p;
  while (isGraph(*p)) ++p;
  if (p - b > width) width = p - b;
  return p;
}

// Copy the next field of src to dest, left-justified
static void sprintLeft(char*& dest, const char*& src, int width) {
  const char *end = dest + width;
  while (isSpace(*src)) ++src;
  while (isGraph(*src)) *dest++ = *src++;
  while (dest < end) *dest++ = ' ';
  *dest++ = ' ';
}

// Copy the next field of src to dest, right-justified
static void sprintRight(char*& dest, const char*& src, int width) {
  while (isSpace(*src)) ++src;
  const char* s = src;
  while (isGraph(*s)) ++s;
  int w = s - src;
  while (w++ < width) *dest++ = ' ';
  while (isGraph(*src)) *dest++ = *src++;
  *dest++ = ' ';
}

static void sprintReplaceRight(char *&dest, const char *&oldIn,
			       const char *newIn, int width) {
  sprintRight(dest, newIn, width);
  oldIn = skipWord(oldIn);
}

// Probability -> phred score in fastq-sanger ASCII representation
static char asciiFromProb(double probRight) {
  double probWrong = 1 - probRight;
  double e = std::max(probWrong, 1e-10);  // avoid overflow errors
  int s = floor(-10 * log10(e));  // phred score, rounded down
  return std::min(s + 33, 126);
}

size_t mafSlice(std::vector<char> &outputText, const UnsplitAlignment &aln,
		unsigned alnBeg, unsigned alnEnd, const double *probs) {
  char begTexts[2][32];
  char lenTexts[2][32];
  int w[6] = {0};
  unsigned numOfLines = 1;  // for an extra "p" line at the end

  int j = 0;
  for (char **i = aln.linesBeg; i < aln.linesEnd; ++i) {
    const char *c = i[0];
    numOfLines += (*c == 's' || *c == 'q' || *c == 'p');
    if (*c == 's') {
      unsigned beg = 0;
      unsigned len = 0;
      c = updateFieldWidth(c, w[0]);
      c = updateFieldWidth(c, w[1]);
      ++c;  // skip over the string terminator
      c = readUint(c, beg);
      c = readUint(c, len);
      c = updateFieldWidth(c, w[4]);
      c = updateFieldWidth(c, w[5]);
      c = skipSpace(c);
      unsigned begPos = seqPosFromAlnPos(alnBeg, c);
      unsigned endPos = seqPosFromAlnPos(alnEnd, c);
      unsigned newBeg = aln.isFlipped() ? beg + len - endPos : beg + begPos;
      unsigned newLen = endPos - begPos;
      w[2] = std::max(w[2], sprintf(begTexts[j], "%u", newBeg));
      w[3] = std::max(w[3], sprintf(lenTexts[j], "%u", newLen));
      ++j;
    }
  }

  unsigned alnLen = alnEnd - alnBeg;
  size_t lineLength = w[0] + w[1] + w[2] + w[3] + w[4] + w[5] + alnLen + 7;
  size_t outputSize = outputText.size();
  outputText.insert(outputText.end(), numOfLines * lineLength, 0);
  char *out = &outputText[outputSize];

  j = 0;
  for (char **i = aln.linesBeg; i < aln.linesEnd; ++i) {
    const char *in = i[0];
    if (*in == 's') {
      sprintLeft(out, in, w[0]);
      sprintLeft(out, in, w[1]);
      ++in;  // skip over the string terminator
      sprintReplaceRight(out, in, begTexts[j], w[2]);
      sprintReplaceRight(out, in, lenTexts[j], w[3]);
      sprintLeft(out, in, w[4]);
      sprintRight(out, in, w[5]);
      in = skipSpace(in);
      memcpy(out, in + alnBeg, alnLen);
      if (aln.isFlipped()) reverseComplement(out, out + alnLen);
      out += alnLen;
      *out++ = '\n';
      ++j;
    } else if (*in == 'q') {
      sprintLeft(out, in, w[0]);
      sprintLeft(out, in, w[1] + w[2] + w[3] + w[4] + w[5] + 4);
      in = skipSpace(in);
      memcpy(out, in + alnBeg, alnLen);
      if (aln.isFlipped()) std::reverse(out, out + alnLen);
      out += alnLen;
      *out++ = '\n';
    } else if (*in == 'p') {
      sprintLeft(out, in, w[0] + w[1] + w[2] + w[3] + w[4] + w[5] + 5);
      const char *beg = aln.isFlipped()
	? rskipSpace(i[1] - 1) - alnEnd
	: skipSpace(in) + alnBeg;
      memcpy(out, beg, alnLen);
      out += alnLen;
      *out++ = '\n';
    }
  }

  const char *in = "p";
  sprintLeft(out, in, w[0] + w[1] + w[2] + w[3] + w[4] + w[5] + 5);
  std::transform(probs, probs + alnLen, out, asciiFromProb);
  if (aln.isFlipped()) std::reverse(out, out + alnLen);
  out += alnLen;
  *out++ = '\n';

  return lineLength;
}

double pLinesToErrorProb(const char *line1, const char *line2) {
  double maxGoodProb = 0;
  const char *i = line1;
  const char *j = line2;
  while (isGraph(*i) && isGraph(*j)) {
    double x = pow(0.1, (*i - 33) * 0.1);  // error probability
    double y = pow(0.1, (*j - 33) * 0.1);  // error probability
    double z = (1 - x) * (1 - y);  // probability that neither is an error
    if (z > maxGoodProb) maxGoodProb = z;
    ++i;
    ++j;
  }
  return 1 - maxGoodProb;  // minimum combined error probability
}

}
