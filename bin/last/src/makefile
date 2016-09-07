CXX = g++
CC  = gcc

CXXFLAGS = -O3 -Wall -Wextra -Wcast-qual -Wswitch-enum -Wundef	\
-Wcast-align -pedantic -g -std=c++11 -pthread -DHAS_CXX_THREADS
# -Wconversion
# -fomit-frame-pointer ?

CFLAGS = -Wall -O2

DBOBJ = Alphabet.o MultiSequence.o CyclicSubsetSeed.o			\
SubsetSuffixArray.o LastdbArguments.o io.o fileMap.o TantanMasker.o	\
ScoreMatrix.o LambdaCalculator.o tantan.o SubsetSuffixArraySort.o	\
MultiSequenceQual.o lastdb.o

ALOBJ = Alphabet.o MultiSequence.o CyclicSubsetSeed.o			\
SubsetSuffixArray.o LastalArguments.o io.o fileMap.o TantanMasker.o	\
ScoreMatrix.o tantan.o DiagonalTable.o SegmentPair.o Alignment.o	\
GappedXdropAligner.o SegmentPairPot.o AlignmentPot.o			\
GeneralizedAffineGapCosts.o Centroid.o LambdaCalculator.o		\
TwoQualityScoreMatrix.o OneQualityScoreMatrix.o QualityPssmMaker.o	\
GeneticCode.o LastEvaluer.o gaplessXdrop.o gaplessPssmXdrop.o		\
gaplessTwoQualityXdrop.o SubsetSuffixArraySearch.o AlignmentWrite.o	\
MultiSequenceQual.o GappedXdropAlignerPssm.o				\
GappedXdropAligner2qual.o GappedXdropAligner3frame.o lastal.o		\
alp/sls_alignment_evaluer.o alp/sls_pvalues.o alp/sls_alp_sim.o		\
alp/sls_alp_regression.o alp/sls_alp_data.o alp/sls_alp.o		\
alp/sls_basic.o alp/njn_localmaxstatmatrix.o alp/njn_localmaxstat.o	\
alp/njn_localmaxstatutil.o alp/njn_dynprogprob.o			\
alp/njn_dynprogprobproto.o alp/njn_dynprogproblim.o alp/njn_ioutil.o	\
alp/njn_random.o alp/sls_falp_alignment_evaluer.o			\
alp/sls_fsa1_pvalues.o alp/sls_fsa1_utils.o alp/sls_fsa1.o		\
alp/sls_fsa1_parameters.o

SPOBJ = Alphabet.o MultiSequence.o fileMap.o split/cbrc_linalg.o	\
split/last-split.o split/cbrc_split_aligner.o split/last-split-main.o	\
split/cbrc_unsplit_alignment.o

PPOBJ = last-pair-probs.o last-pair-probs-main.o io.o

MBOBJ = last-merge-batches.o

ALL = lastdb lastal last-split last-merge-batches last-pair-probs

all: $(ALL)

lastdb: $(DBOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(DBOBJ)

lastal: $(ALOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(ALOBJ)

last-split: $(SPOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SPOBJ)

last-pair-probs: $(PPOBJ)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $(PPOBJ)

last-merge-batches: $(MBOBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(MBOBJ)

.SUFFIXES:
.SUFFIXES: .o .c .cc .cpp

.c.o:
	$(CC) $(CPPFLAGS) $(CFLAGS) -c -o $@ $<

.cc.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I. -c -o $@ $<

.cpp.o:
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -f $(ALL) *.o */*.o

CyclicSubsetSeedData.hh: ../data/*.seed
	../build/seed-inc.sh ../data/*.seed > $@

ScoreMatrixData.hh: ../data/*.mat
	../build/mat-inc.sh ../data/*.mat > $@

VERSION = \"`hg id -n`\"
VERSION2 = \"`grep latesttagdistance ../.hg_archival.txt | cut -d' ' -f2`\"
VERSION3 = \"UNKNOWN\"

version.hh: FORCE
	if test -e ../.hg ; \
	then echo $(VERSION) | cmp -s $@ - || echo $(VERSION) > $@ ; \
	elif test -e ../.hg_archival.txt ; \
	then echo $(VERSION2) | cmp -s $@ - || echo $(VERSION2) > $@ ; \
	else test -e $@ || echo $(VERSION3) > $@ ; \
	fi

FORCE:

depend:
	sed '/[m][v]/q' makefile > m
	$(CXX) -MM *.cc >> m
	$(CC) -MM *.c >> m
	$(CXX) -MM alp/*.cpp | sed 's|.*:|alp/&|' >> m
	$(CXX) -MM -I. split/*.cc | sed 's|.*:|split/&|' >> m
	mv m makefile
Alignment.o: Alignment.cc Alignment.hh ScoreMatrixRow.hh SegmentPair.hh \
 Alphabet.hh Centroid.hh GappedXdropAligner.hh \
 GeneralizedAffineGapCosts.hh OneQualityScoreMatrix.hh GeneticCode.hh \
 TwoQualityScoreMatrix.hh
AlignmentPot.o: AlignmentPot.cc AlignmentPot.hh Alignment.hh \
 ScoreMatrixRow.hh SegmentPair.hh
AlignmentWrite.o: AlignmentWrite.cc Alignment.hh ScoreMatrixRow.hh \
 SegmentPair.hh GeneticCode.hh LastEvaluer.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 MultiSequence.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 Alphabet.hh
Alphabet.o: Alphabet.cc Alphabet.hh
Centroid.o: Centroid.cc Centroid.hh GappedXdropAligner.hh \
 ScoreMatrixRow.hh GeneralizedAffineGapCosts.hh SegmentPair.hh \
 OneQualityScoreMatrix.hh GappedXdropAlignerInl.hh
CyclicSubsetSeed.o: CyclicSubsetSeed.cc CyclicSubsetSeed.hh \
 CyclicSubsetSeedData.hh io.hh stringify.hh
DiagonalTable.o: DiagonalTable.cc DiagonalTable.hh
GappedXdropAligner.o: GappedXdropAligner.cc GappedXdropAligner.hh \
 ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAligner2qual.o: GappedXdropAligner2qual.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh \
 TwoQualityScoreMatrix.hh
GappedXdropAligner3frame.o: GappedXdropAligner3frame.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAligner3framePssm.o: GappedXdropAligner3framePssm.cc \
 GappedXdropAligner.hh ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GappedXdropAlignerPssm.o: GappedXdropAlignerPssm.cc GappedXdropAligner.hh \
 ScoreMatrixRow.hh GappedXdropAlignerInl.hh
GeneralizedAffineGapCosts.o: GeneralizedAffineGapCosts.cc \
 GeneralizedAffineGapCosts.hh
GeneticCode.o: GeneticCode.cc GeneticCode.hh Alphabet.hh
LambdaCalculator.o: LambdaCalculator.cc LambdaCalculator.hh
LastEvaluer.o: LastEvaluer.cc LastEvaluer.hh ScoreMatrixRow.hh \
 alp/sls_alignment_evaluer.hpp alp/sls_pvalues.hpp alp/sls_basic.hpp \
 alp/sls_falp_alignment_evaluer.hpp alp/sls_fsa1_pvalues.hpp \
 GeneticCode.hh
LastalArguments.o: LastalArguments.cc LastalArguments.hh \
 SequenceFormat.hh stringify.hh version.hh
LastdbArguments.o: LastdbArguments.cc LastdbArguments.hh \
 SequenceFormat.hh stringify.hh version.hh
MultiSequence.o: MultiSequence.cc MultiSequence.hh ScoreMatrixRow.hh \
 VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh io.hh
MultiSequenceQual.o: MultiSequenceQual.cc MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
OneQualityScoreMatrix.o: OneQualityScoreMatrix.cc \
 OneQualityScoreMatrix.hh ScoreMatrixRow.hh qualityScoreUtil.hh \
 stringify.hh
QualityPssmMaker.o: QualityPssmMaker.cc QualityPssmMaker.hh \
 ScoreMatrixRow.hh qualityScoreUtil.hh stringify.hh
ScoreMatrix.o: ScoreMatrix.cc ScoreMatrix.hh ScoreMatrixData.hh io.hh
SegmentPair.o: SegmentPair.cc SegmentPair.hh
SegmentPairPot.o: SegmentPairPot.cc SegmentPairPot.hh SegmentPair.hh
SubsetSuffixArray.o: SubsetSuffixArray.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 io.hh
SubsetSuffixArraySearch.o: SubsetSuffixArraySearch.cc \
 SubsetSuffixArray.hh CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh
SubsetSuffixArraySort.o: SubsetSuffixArraySort.cc SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
TantanMasker.o: TantanMasker.cc TantanMasker.hh ScoreMatrixRow.hh \
 tantan.hh ScoreMatrix.hh LambdaCalculator.hh
TwoQualityScoreMatrix.o: TwoQualityScoreMatrix.cc \
 TwoQualityScoreMatrix.hh ScoreMatrixRow.hh qualityScoreUtil.hh \
 stringify.hh
fileMap.o: fileMap.cc fileMap.hh stringify.hh
gaplessPssmXdrop.o: gaplessPssmXdrop.cc gaplessPssmXdrop.hh \
 ScoreMatrixRow.hh
gaplessTwoQualityXdrop.o: gaplessTwoQualityXdrop.cc \
 gaplessTwoQualityXdrop.hh TwoQualityScoreMatrix.hh ScoreMatrixRow.hh
gaplessXdrop.o: gaplessXdrop.cc gaplessXdrop.hh ScoreMatrixRow.hh
io.o: io.cc io.hh
last-pair-probs-main.o: last-pair-probs-main.cc last-pair-probs.hh \
 stringify.hh version.hh
last-pair-probs.o: last-pair-probs.cc last-pair-probs.hh io.hh \
 stringify.hh
lastal.o: lastal.cc LastalArguments.hh SequenceFormat.hh \
 QualityPssmMaker.hh ScoreMatrixRow.hh OneQualityScoreMatrix.hh \
 TwoQualityScoreMatrix.hh qualityScoreUtil.hh stringify.hh \
 LambdaCalculator.hh LastEvaluer.hh alp/sls_alignment_evaluer.hpp \
 alp/sls_pvalues.hpp alp/sls_basic.hpp alp/sls_falp_alignment_evaluer.hpp \
 alp/sls_fsa1_pvalues.hpp GeneticCode.hh SubsetSuffixArray.hh \
 CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh fileMap.hh Centroid.hh \
 GappedXdropAligner.hh GeneralizedAffineGapCosts.hh SegmentPair.hh \
 AlignmentPot.hh Alignment.hh SegmentPairPot.hh ScoreMatrix.hh \
 Alphabet.hh MultiSequence.hh TantanMasker.hh tantan.hh DiagonalTable.hh \
 gaplessXdrop.hh gaplessPssmXdrop.hh gaplessTwoQualityXdrop.hh io.hh \
 version.hh
lastdb.o: lastdb.cc LastdbArguments.hh SequenceFormat.hh \
 SubsetSuffixArray.hh CyclicSubsetSeed.hh VectorOrMmap.hh Mmap.hh \
 fileMap.hh stringify.hh Alphabet.hh MultiSequence.hh ScoreMatrixRow.hh \
 TantanMasker.hh tantan.hh io.hh qualityScoreUtil.hh version.hh
tantan.o: tantan.cc tantan.hh
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
split/cbrc_linalg.o: split/cbrc_linalg.cc split/cbrc_linalg.hh
split/cbrc_split_aligner.o: split/cbrc_split_aligner.cc \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh \
 split/cbrc_linalg.hh
split/cbrc_unsplit_alignment.o: split/cbrc_unsplit_alignment.cc \
 split/cbrc_unsplit_alignment.hh
split/last-split-main.o: split/last-split-main.cc split/last-split.hh \
 stringify.hh version.hh
split/last-split.o: split/last-split.cc split/last-split.hh \
 split/cbrc_split_aligner.hh split/cbrc_unsplit_alignment.hh \
 split/cbrc_int_exponentiator.hh Alphabet.hh MultiSequence.hh \
 ScoreMatrixRow.hh VectorOrMmap.hh Mmap.hh fileMap.hh stringify.hh
