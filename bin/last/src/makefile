CXXFLAGS = -msse4 -O3 -Wall -g -std=c++11 -pthread
# -fomit-frame-pointer ?

# The number of symbol types that we can keep in sequences (e.g. 20
# amino acids, plus ambiguous ones, in upper & lower case, plus one
# delimiter):
ALPHABET_CAPACITY = 66
CPPF = -DALPHABET_CAPACITY=$(ALPHABET_CAPACITY) -DHAS_CXX_THREADS $(CPPFLAGS)

CFLAGS = -Wall -O2

alpObj = alp/sls_alignment_evaluer.o alp/sls_pvalues.o		\
alp/sls_alp_sim.o alp/sls_alp_regression.o alp/sls_alp_data.o	\
alp/sls_alp.o alp/sls_basic.o alp/njn_localmaxstatmatrix.o	\
alp/njn_localmaxstat.o alp/njn_localmaxstatutil.o		\
alp/njn_dynprogprob.o alp/njn_dynprogprobproto.o		\
alp/njn_dynprogproblim.o alp/njn_ioutil.o alp/njn_random.o	\
alp/sls_falp_alignment_evaluer.o alp/sls_fsa1_pvalues.o		\
alp/sls_fsa1_utils.o alp/sls_fsa1.o alp/sls_fsa1_parameters.o

indexObj0 = Alphabet.o CyclicSubsetSeed.o LambdaCalculator.o	\
ScoreMatrix.o SubsetMinimizerFinder.o TantanMasker.o		\
dna_words_finder.o fileMap.o cbrc_linalg.o			\
mcf_substitution_matrix_stats.o tantan.o LastdbArguments.o

indexObj4 = MultiSequence.o MultiSequenceQual.o SubsetSuffixArray.o	\
SubsetSuffixArraySort.o lastdb.o

alignObj0 = Alphabet.o CyclicSubsetSeed.o LambdaCalculator.o		\
ScoreMatrix.o SubsetMinimizerFinder.o TantanMasker.o			\
dna_words_finder.o fileMap.o tantan.o LastalArguments.o			\
GappedXdropAligner.o GappedXdropAlignerDna.o GappedXdropAlignerPssm.o	\
GappedXdropAligner2qual.o GappedXdropAligner3frame.o			\
GappedXdropAlignerFrame.o mcf_aligment_path_adder.o			\
mcf_frameshift_xdrop_aligner.o mcf_gap_costs.o GeneticCode.o		\
GreedyXdropAligner.o LastEvaluer.o OneQualityScoreMatrix.o		\
QualityPssmMaker.o TwoQualityScoreMatrix.o gaplessXdrop.o		\
gaplessPssmXdrop.o gaplessTwoQualityXdrop.o cbrc_linalg.o		\
mcf_substitution_matrix_stats.o split/cbrc_unsplit_alignment.o		\
split/last_split_options.o $(alpObj)

alignObj4 = Alignment.o AlignmentPot.o AlignmentWrite.o Centroid.o	\
DiagonalTable.o MultiSequence.o MultiSequenceQual.o SegmentPair.o	\
SegmentPairPot.o SubsetSuffixArray.o SubsetSuffixArraySearch.o		\
lastal.o split/cbrc_split_aligner.o split/mcf_last_splitter.o

splitObj0 = Alphabet.o LambdaCalculator.o fileMap.o cbrc_linalg.o	\
mcf_substitution_matrix_stats.o split/cbrc_unsplit_alignment.o		\
split/last_split_options.o split/last-split-main.o

splitObj4 = MultiSequence.o split/cbrc_split_aligner.o	\
split/mcf_last_splitter.o split/last-split.o

PPOBJ = last-pair-probs.o last-pair-probs-main.o

MBOBJ = last-merge-batches.o

LAST4 = ../bin/lastdb ../bin/lastal ../bin/last-split
LAST5 = ../bin/lastdb5 ../bin/lastal5 ../bin/last-split5
LAST8 = ../bin/lastdb8 ../bin/lastal8 ../bin/last-split8
ALL = $(LAST4) $(LAST5) ../bin/last-merge-batches ../bin/last-pair-probs

indexObj5 = $(indexObj4:.o=.o5)
alignObj5 = $(alignObj4:.o=.o5)
splitObj5 = $(splitObj4:.o=.o5)

indexObj8 = $(indexObj4:.o=.o8)
alignObj8 = $(alignObj4:.o=.o8)
splitObj8 = $(splitObj4:.o=.o8)

all: $(ALL)
last8: $(LAST8)

indexAllObj4 = $(indexObj0) $(indexObj4)
../bin/lastdb: $(indexAllObj4)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(indexAllObj4) -lz

indexAllObj5 = $(indexObj0) $(indexObj5)
../bin/lastdb5: $(indexAllObj5)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(indexAllObj5) -lz

indexAllObj8 = $(indexObj0) $(indexObj8)
../bin/lastdb8: $(indexAllObj8)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(indexAllObj8) -lz

alignAllObj4 = $(alignObj0) $(alignObj4)
../bin/lastal: $(alignAllObj4)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(alignAllObj4) -lz

alignAllObj5 = $(alignObj0) $(alignObj5)
../bin/lastal5: $(alignAllObj5)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(alignAllObj5) -lz

alignAllObj8 = $(alignObj0) $(alignObj8)
../bin/lastal8: $(alignAllObj8)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(alignAllObj8) -lz

splitAllObj4 = $(splitObj0) $(splitObj4)
../bin/last-split: $(splitAllObj4)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(splitAllObj4)

splitAllObj5 = $(splitObj0) $(splitObj5)
../bin/last-split5: $(splitAllObj5)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(splitAllObj5)

splitAllObj8 = $(splitObj0) $(splitObj8)
../bin/last-split8: $(splitAllObj8)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(splitAllObj8)

../bin/last-pair-probs: $(PPOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(PPOBJ) -lz

../bin/last-merge-batches: $(MBOBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(MBOBJ)

.SUFFIXES:
.SUFFIXES: .o .c .cc .cpp .o5 .o8

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

.cc.o:
	$(CXX) -DLAST_POS_BYTES=4 $(CPPF) $(CXXFLAGS) -I. -c -o $@ $<

.cc.o5:
	$(CXX) -DLAST_POS_BYTES=5 $(CPPF) $(CXXFLAGS) -I. -c -o $@ $<

.cc.o8:
	$(CXX) -DLAST_POS_BYTES=8 $(CPPF) $(CXXFLAGS) -I. -c -o $@ $<

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(ALL) $(LAST8) *.o* */*.o*

CyclicSubsetSeedData.hh: ../data/*.seed
	../build/seed-inc.sh ../data/*.seed > $@

GeneticCodeData.hh: ../data/gc.prt
	../build/gc-inc.sh ../data/gc.prt > $@

ScoreMatrixData.hh: ../data/*.mat
	../build/mat-inc.sh ../data/*.mat > $@

VERSION1 = git describe --dirty
VERSION2 = echo ' (HEAD -> main, tag: 1418) ' | sed -e 's/.*tag: *//' -e 's/[,) ].*//'

VERSION = \"`test -e ../.git && $(VERSION1) || $(VERSION2)`\"

version.hh: FORCE
	echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@

FORCE:

fix8 = sed 's/.*\.o/& &5 &8/'

depend:
	sed '/[m][v]/q' makefile > m
	$(CXX) -MM -I. -std=c++11 *.cc | $(fix8) >> m
	$(CC) -MM *.c >> m
	$(CXX) -MM alp/*.cpp | sed 's|.*:|alp/&|' >> m
	$(CXX) -MM -I. split/*.cc | sed 's|.*:|split/&|' | $(fix8) >> m
	mv m makefile
Alignment.o Alignment.o5 Alignment.o8: Alignment.cc Alignment.hh Centroid.hh GappedXdropAligner.hh \
 mcf_contiguous_queue.hh mcf_reverse_queue.hh mcf_gap_costs.hh \
 mcf_simd.hh ScoreMatrixRow.hh SegmentPair.hh OneQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh GreedyXdropAligner.hh \
 mcf_frameshift_xdrop_aligner.hh Alphabet.hh GeneticCode.hh \
 TwoQualityScoreMatrix.hh
AlignmentPot.o AlignmentPot.o5 AlignmentPot.o8: AlignmentPot.cc AlignmentPot.hh Alignment.hh Centroid.hh \
 GappedXdropAligner.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh SegmentPair.hh \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 GreedyXdropAligner.hh mcf_frameshift_xdrop_aligner.hh
AlignmentWrite.o AlignmentWrite.o5 AlignmentWrite.o8: AlignmentWrite.cc Alignment.hh Centroid.hh \
 GappedXdropAligner.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh SegmentPair.hh \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 GreedyXdropAligner.hh mcf_frameshift_xdrop_aligner.hh GeneticCode.hh \
 LastEvaluer.hh alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp \
 alp/sls_basic.hpp MultiSequence.hh VectorOrMmap.hh Mmap.hh fileMap.hh \
 stringify.hh Alphabet.hh
Alphabet.o Alphabet.o5 Alphabet.o8: Alphabet.cc Alphabet.hh
cbrc_linalg.o cbrc_linalg.o5 cbrc_linalg.o8: cbrc_linalg.cc cbrc_linalg.hh
Centroid.o Centroid.o5 Centroid.o8: Centroid.cc Centroid.hh GappedXdropAligner.hh \
 mcf_contiguous_queue.hh mcf_reverse_queue.hh mcf_gap_costs.hh \
 mcf_simd.hh ScoreMatrixRow.hh SegmentPair.hh OneQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh GappedXdropAlignerInl.hh
CyclicSubsetSeed.o CyclicSubsetSeed.o5 CyclicSubsetSeed.o8: CyclicSubsetSeed.cc CyclicSubsetSeed.hh \
 CyclicSubsetSeedData.hh zio.hh mcf_zstream.hh stringify.hh
DiagonalTable.o DiagonalTable.o5 DiagonalTable.o8: DiagonalTable.cc DiagonalTable.hh
dna_words_finder.o dna_words_finder.o5 dna_words_finder.o8: dna_words_finder.cc dna_words_finder.hh
fileMap.o fileMap.o5 fileMap.o8: fileMap.cc fileMap.hh stringify.hh
gaplessPssmXdrop.o gaplessPssmXdrop.o5 gaplessPssmXdrop.o8: gaplessPssmXdrop.cc gaplessPssmXdrop.hh \
 ScoreMatrixRow.hh
gaplessTwoQualityXdrop.o gaplessTwoQualityXdrop.o5 gaplessTwoQualityXdrop.o8: gaplessTwoQualityXdrop.cc \
 gaplessTwoQualityXdrop.hh TwoQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh ScoreMatrixRow.hh
gaplessXdrop.o gaplessXdrop.o5 gaplessXdrop.o8: gaplessXdrop.cc gaplessXdrop.hh ScoreMatrixRow.hh
GappedXdropAligner2qual.o GappedXdropAligner2qual.o5 GappedXdropAligner2qual.o8: GappedXdropAligner2qual.cc \
 GappedXdropAligner.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh \
 TwoQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh
GappedXdropAligner3frame.o GappedXdropAligner3frame.o5 GappedXdropAligner3frame.o8: GappedXdropAligner3frame.cc \
 GappedXdropAligner.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAligner3framePssm.o GappedXdropAligner3framePssm.o5 GappedXdropAligner3framePssm.o8: GappedXdropAligner3framePssm.cc \
 GappedXdropAligner.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAligner.o GappedXdropAligner.o5 GappedXdropAligner.o8: GappedXdropAligner.cc GappedXdropAligner.hh \
 mcf_contiguous_queue.hh mcf_reverse_queue.hh mcf_gap_costs.hh \
 mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAlignerDna.o GappedXdropAlignerDna.o5 GappedXdropAlignerDna.o8: GappedXdropAlignerDna.cc GappedXdropAligner.hh \
 mcf_contiguous_queue.hh mcf_reverse_queue.hh mcf_gap_costs.hh \
 mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAlignerFrame.o GappedXdropAlignerFrame.o5 GappedXdropAlignerFrame.o8: GappedXdropAlignerFrame.cc \
 GappedXdropAligner.hh mcf_contiguous_queue.hh mcf_reverse_queue.hh \
 mcf_gap_costs.hh mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAlignerPssm.o GappedXdropAlignerPssm.o5 GappedXdropAlignerPssm.o8: GappedXdropAlignerPssm.cc GappedXdropAligner.hh \
 mcf_contiguous_queue.hh mcf_reverse_queue.hh mcf_gap_costs.hh \
 mcf_simd.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GeneticCode.o GeneticCode.o5 GeneticCode.o8: GeneticCode.cc GeneticCode.hh GeneticCodeData.hh \
 Alphabet.hh zio.hh mcf_zstream.hh
GreedyXdropAligner.o GreedyXdropAligner.o5 GreedyXdropAligner.o8: GreedyXdropAligner.cc GreedyXdropAligner.hh \
 ScoreMatrixRow.hh
LambdaCalculator.o LambdaCalculator.o5 LambdaCalculator.o8: LambdaCalculator.cc LambdaCalculator.hh \
 cbrc_linalg.hh
LastalArguments.o LastalArguments.o5 LastalArguments.o8: LastalArguments.cc LastalArguments.hh \
 SequenceFormat.hh split/last_split_options.hh stringify.hh getoptUtil.hh \
 version.hh
lastal.o lastal.o5 lastal.o8: lastal.cc last.hh Alphabet.hh CyclicSubsetSeed.hh \
 MultiSequence.hh ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh \
 stringify.hh SequenceFormat.hh SubsetSuffixArray.hh dna_words_finder.hh \
 qualityScoreUtil.hh LastalArguments.hh split/last_split_options.hh \
 QualityPssmMaker.hh OneQualityScoreMatrix.hh \
 mcf_substitution_matrix_stats.hh TwoQualityScoreMatrix.hh LastEvaluer.hh \
 mcf_frameshift_xdrop_aligner.hh mcf_gap_costs.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 GeneticCode.hh SubsetMinimizerFinder.hh AlignmentPot.hh Alignment.hh \
 Centroid.hh GappedXdropAligner.hh mcf_contiguous_queue.hh \
 mcf_reverse_queue.hh mcf_simd.hh SegmentPair.hh GreedyXdropAligner.hh \
 SegmentPairPot.hh ScoreMatrix.hh TantanMasker.hh tantan.hh \
 DiagonalTable.hh gaplessXdrop.hh gaplessPssmXdrop.hh \
 gaplessTwoQualityXdrop.hh zio.hh mcf_zstream.hh threadUtil.hh \
 split/mcf_last_splitter.hh split/cbrc_split_aligner.hh \
 split/cbrc_unsplit_alignment.hh split/cbrc_int_exponentiator.hh \
 Alphabet.hh MultiSequence.hh split/last_split_options.hh version.hh
LastdbArguments.o LastdbArguments.o5 LastdbArguments.o8: LastdbArguments.cc LastdbArguments.hh \
 SequenceFormat.hh stringify.hh getoptUtil.hh version.hh
lastdb.o lastdb.o5 lastdb.o8: lastdb.cc last.hh Alphabet.hh CyclicSubsetSeed.hh \
 MultiSequence.hh ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh \
 stringify.hh SequenceFormat.hh SubsetSuffixArray.hh dna_words_finder.hh \
 qualityScoreUtil.hh LastdbArguments.hh TantanMasker.hh tantan.hh zio.hh \
 mcf_zstream.hh threadUtil.hh version.hh
LastEvaluer.o LastEvaluer.o5 LastEvaluer.o8: LastEvaluer.cc LastEvaluer.hh ScoreMatrixRow.hh \
 mcf_frameshift_xdrop_aligner.hh mcf_gap_costs.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 GeneticCode.hh mcf_aligment_path_adder.hh \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp
last-pair-probs.o last-pair-probs.o5 last-pair-probs.o8: last-pair-probs.cc last-pair-probs.hh zio.hh \
 mcf_zstream.hh stringify.hh
last-pair-probs-main.o last-pair-probs-main.o5 last-pair-probs-main.o8: last-pair-probs-main.cc last-pair-probs.hh \
 stringify.hh version.hh
mcf_aligment_path_adder.o mcf_aligment_path_adder.o5 mcf_aligment_path_adder.o8: mcf_aligment_path_adder.cc \
 mcf_aligment_path_adder.hh
mcf_frameshift_xdrop_aligner.o mcf_frameshift_xdrop_aligner.o5 mcf_frameshift_xdrop_aligner.o8: mcf_frameshift_xdrop_aligner.cc \
 mcf_frameshift_xdrop_aligner.hh mcf_gap_costs.hh
mcf_gap_costs.o mcf_gap_costs.o5 mcf_gap_costs.o8: mcf_gap_costs.cc mcf_gap_costs.hh
mcf_substitution_matrix_stats.o mcf_substitution_matrix_stats.o5 mcf_substitution_matrix_stats.o8: mcf_substitution_matrix_stats.cc \
 mcf_substitution_matrix_stats.hh LambdaCalculator.hh cbrc_linalg.hh
MultiSequence.o MultiSequence.o5 MultiSequence.o8: MultiSequence.cc MultiSequence.hh ScoreMatrixRow.hh \
 VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh io.hh
MultiSequenceQual.o MultiSequenceQual.o5 MultiSequenceQual.o8: MultiSequenceQual.cc MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
OneQualityScoreMatrix.o OneQualityScoreMatrix.o5 OneQualityScoreMatrix.o8: OneQualityScoreMatrix.cc \
 OneQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
QualityPssmMaker.o QualityPssmMaker.o5 QualityPssmMaker.o8: QualityPssmMaker.cc QualityPssmMaker.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
ScoreMatrix.o ScoreMatrix.o5 ScoreMatrix.o8: ScoreMatrix.cc ScoreMatrix.hh ScoreMatrixRow.hh \
 ScoreMatrixData.hh qualityScoreUtil.hh stringify.hh zio.hh \
 mcf_zstream.hh
SegmentPair.o SegmentPair.o5 SegmentPair.o8: SegmentPair.cc SegmentPair.hh
SegmentPairPot.o SegmentPairPot.o5 SegmentPairPot.o8: SegmentPairPot.cc SegmentPairPot.hh SegmentPair.hh
SubsetMinimizerFinder.o SubsetMinimizerFinder.o5 SubsetMinimizerFinder.o8: SubsetMinimizerFinder.cc \
 SubsetMinimizerFinder.hh CyclicSubsetSeed.hh
SubsetSuffixArray.o SubsetSuffixArray.o5 SubsetSuffixArray.o8: SubsetSuffixArray.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh dna_words_finder.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh SubsetMinimizerFinder.hh io.hh
SubsetSuffixArraySearch.o SubsetSuffixArraySearch.o5 SubsetSuffixArraySearch.o8: SubsetSuffixArraySearch.cc \
 SubsetSuffixArray.hh CyclicSubsetSeed.hh dna_words_finder.hh \
 VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
SubsetSuffixArraySort.o SubsetSuffixArraySort.o5 SubsetSuffixArraySort.o8: SubsetSuffixArraySort.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh dna_words_finder.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh
tantan.o tantan.o5 tantan.o8: tantan.cc tantan.hh mcf_simd.hh
TantanMasker.o TantanMasker.o5 TantanMasker.o8: TantanMasker.cc TantanMasker.hh ScoreMatrixRow.hh \
 tantan.hh ScoreMatrix.hh mcf_substitution_matrix_stats.hh
TwoQualityScoreMatrix.o TwoQualityScoreMatrix.o5 TwoQualityScoreMatrix.o8: TwoQualityScoreMatrix.cc \
 TwoQualityScoreMatrix.hh mcf_substitution_matrix_stats.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
last-merge-batches.o: last-merge-batches.c version.hh
alp/njn_dynprogprob.o: alp/njn_dynprogprob.cpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp
alp/njn_dynprogproblim.o: alp/njn_dynprogproblim.cpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp
alp/njn_dynprogprobproto.o: alp/njn_dynprogprobproto.cpp \
 alp/njn_dynprogprobproto.hpp
alp/njn_ioutil.o: alp/njn_ioutil.cpp alp/njn_ioutil.hpp
alp/njn_localmaxstat.o: alp/njn_localmaxstat.cpp alp/sls_basic.hpp \
 alp/njn_localmaxstat.hpp alp/njn_memutil.hpp alp/njn_ioutil.hpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_function.hpp alp/njn_doubletype.hpp \
 alp/njn_integer.hpp alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp \
 alp/njn_approx.hpp alp/njn_vector.hpp
alp/njn_localmaxstatmatrix.o: alp/njn_localmaxstatmatrix.cpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/njn_memutil.hpp
alp/njn_localmaxstatutil.o: alp/njn_localmaxstatutil.cpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/njn_dynprogproblim.hpp alp/njn_dynprogprob.hpp \
 alp/njn_dynprogprobproto.hpp alp/njn_integer.hpp alp/njn_memutil.hpp \
 alp/njn_root.hpp alp/njn_function.hpp alp/sls_basic.hpp
alp/njn_random.o: alp/njn_random.cpp alp/njn_random.hpp
alp/sls_alignment_evaluer.o: alp/sls_alignment_evaluer.cpp \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_alp.hpp alp/sls_alp_data.hpp alp/sls_alp_regression.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_sim.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp
alp/sls_alp.o: alp/sls_alp.cpp alp/sls_alp.hpp alp/sls_alp_data.hpp \
 alp/sls_basic.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp
alp/sls_alp_data.o: alp/sls_alp_data.cpp alp/sls_alp_data.hpp \
 alp/sls_basic.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp
alp/sls_alp_regression.o: alp/sls_alp_regression.cpp \
 alp/sls_alp_regression.hpp alp/sls_basic.hpp
alp/sls_alp_sim.o: alp/sls_alp_sim.cpp alp/sls_alp_sim.hpp \
 alp/sls_alp_data.hpp alp/sls_basic.hpp alp/sls_alp_regression.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp.hpp
alp/sls_basic.o: alp/sls_basic.cpp alp/sls_basic.hpp
alp/sls_falp_alignment_evaluer.o: alp/sls_falp_alignment_evaluer.cpp \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_parameters.hpp alp/sls_fsa1_utils.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_regression.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp \
 alp/sls_fsa1.hpp
alp/sls_fsa1.o: alp/sls_fsa1.cpp alp/sls_fsa1.hpp alp/sls_alp_regression.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_utils.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_fsa1_parameters.hpp alp/sls_fsa1_pvalues.hpp \
 alp/njn_localmaxstatmatrix.hpp alp/njn_localmaxstat.hpp \
 alp/njn_localmaxstatutil.hpp alp/njn_matrix.hpp alp/njn_approx.hpp \
 alp/njn_doubletype.hpp alp/njn_ioutil.hpp alp/njn_vector.hpp
alp/sls_fsa1_parameters.o: alp/sls_fsa1_parameters.cpp \
 alp/sls_fsa1_parameters.hpp alp/sls_fsa1_utils.hpp alp/sls_basic.hpp \
 alp/njn_random.hpp alp/njn_uniform.hpp alp/sls_alp_regression.hpp \
 alp/sls_fsa1_pvalues.hpp
alp/sls_fsa1_pvalues.o: alp/sls_fsa1_pvalues.cpp alp/sls_fsa1_pvalues.hpp \
 alp/sls_basic.hpp alp/sls_fsa1_utils.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_normal_distr_array.hpp
alp/sls_fsa1_utils.o: alp/sls_fsa1_utils.cpp alp/sls_fsa1_utils.hpp \
 alp/sls_basic.hpp alp/njn_random.hpp alp/njn_uniform.hpp
alp/sls_pvalues.o: alp/sls_pvalues.cpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_alp_data.hpp alp/sls_alp_regression.hpp alp/njn_random.hpp \
 alp/njn_uniform.hpp alp/sls_normal_distr_array.hpp
split/cbrc_split_aligner.o split/cbrc_split_aligner.o5 split/cbrc_split_aligner.o8: split/cbrc_split_aligner.cc \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 mcf_substitution_matrix_stats.hh
split/cbrc_unsplit_alignment.o split/cbrc_unsplit_alignment.o5 split/cbrc_unsplit_alignment.o8: split/cbrc_unsplit_alignment.cc \
 split/cbrc_unsplit_alignment.hh
split/last-split.o split/last-split.o5 split/last-split.o8: split/last-split.cc split/last-split.hh \
 split/last_split_options.hh split/mcf_last_splitter.hh \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
split/last-split-main.o split/last-split-main.o5 split/last-split-main.o8: split/last-split-main.cc split/last-split.hh \
 split/last_split_options.hh stringify.hh version.hh
split/last_split_options.o split/last_split_options.o5 split/last_split_options.o8: split/last_split_options.cc \
 split/last_split_options.hh
split/mcf_last_splitter.o split/mcf_last_splitter.o5 split/mcf_last_splitter.o8: split/mcf_last_splitter.cc \
 split/mcf_last_splitter.hh split/cbrc_split_aligner.hh \
 split/cbrc_unsplit_alignment.hh split/cbrc_int_exponentiator.hh \
 Alphabet.hh MultiSequence.hh ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh split/last_split_options.hh
