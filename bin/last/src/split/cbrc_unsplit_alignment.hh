// Copyright 2012 Risa Kawaguchi
// Copyright 2013, 2014 Martin C. Frith

//
// 1. UnsplitAlignment -> container of maf file information
//

#ifndef CBRC_UNSPLIT_ALIGNMENT_HH
#define CBRC_UNSPLIT_ALIGNMENT_HH

#include <stddef.h>

#include <vector>

namespace cbrc {

class UnsplitAlignment {
public:
    char **linesBeg;
    char **linesEnd;
    const char *qname;
    unsigned qstart;
    unsigned qend;
    char qstrand;
    unsigned rstart;
    unsigned rend;
    const char *rname;
    const char *ralign;
    const char *qalign;
    const char *qQual;
    UnsplitAlignment(){}
    UnsplitAlignment(char **linesBegIn, char **linesEndIn, bool isTopSeqQuery)
      : linesBeg(linesBegIn), linesEnd(linesEndIn) { init(isTopSeqQuery); }
    void init(bool isTopSeqQuery);
    bool isForwardStrand() const { return qstrand < 2; }
    bool isFlipped() const { return qstrand % 2; }
};

// Appends maf "s", "q", and "p" lines to outputText.
// Appends an extra "p" line for "probs".
// Returns the line length (including a newline).
size_t mafSlice(std::vector<char> &outputText, const UnsplitAlignment &aln,
		unsigned alnBeg, unsigned alnEnd, const double *probs);

void mafSliceBeg(const char* rAln, const char* qAln,
		 unsigned qBeg, unsigned& qSliceBeg, unsigned& alnBeg);

void mafSliceEnd(const char* rAln, const char* qAln,
		 unsigned qEnd, unsigned& qSliceEnd, unsigned& alnEnd);

double pLinesToErrorProb(const char *line1, const char *line2);

}

#endif
