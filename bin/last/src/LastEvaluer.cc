// Copyright 2015 Martin C. Frith

#include "LastEvaluer.hh"

#include "GeneticCode.hh"

#include <algorithm>
#include <cstring>

#define COUNTOF(a) (sizeof (a) / sizeof *(a))

namespace cbrc {

// These lookup tables are not necessary: they are intended to make
// startup faster for typical parameters.

struct EvalueParametersByName {
  const char *matrixName;
  int gapOpen;
  int gapEpen;
  Sls::AlignmentEvaluerParameters parameters;
};

struct EvalueParametersByScore {
  int matchScore;
  int mismatchCost;
  int gapOpen;
  int gapEpen;
  Sls::AlignmentEvaluerParameters parameters;
};

struct FrameshiftEvalueParameters {
  const char *matrixName;
  int gapOpen;
  int gapEpen;
  int frameshiftCost;
  Sls::AlignmentEvaluerParameters parameters;
};

const EvalueParametersByName proteinParameters[] = {
  {"BL62", 11, 1, {0.27553858462416075, 0.052108666891553253,
		   1.5606171794974455, -19.732860493752256,
		   1.5606171794974455, -19.732860493752256,
		   29.210151638988194, -598.91453155494662,
		   29.210151638988194, -598.91453155494662,
		   28.51287191489801, -582.17981817678219}},

  {"BL62", 11, 2, {0.30593548727296388, 0.10240374069870158,
		   0.99050553069382763, -6.5543626660042129,
		   0.99050553069382763, -6.5543626660042129,
		   9.286796977904153, -130.81685466300709,
		   9.286796977904153, -130.81685466300709,
		   9.2266144189146271, -129.25210812927941}},

  {"BL80", 11, 1, {0.32849624141749045, 0.1163738959652578,
		   0.75491225062543932, -6.0768179898445682,
		   0.75491225062543932, -6.0768179898445682,
		   5.5209318716230964, -92.3656298469962,
		   5.5209318716230964, -92.3656298469962,
		   5.2639936937083576, -86.199113577042453}},

  {"BL80", 11, 2, {0.34392164286655347, 0.16475434978940937,
		   0.58112958853087737, -2.0648702745396719,
		   0.58112958853087737, -2.0648702745396719,
		   2.187240038686872, -13.386778011237375,
		   2.187240038686872, -13.386778011237375,
		   2.1552267345409919, -12.554432103444491}},

  {"MIQS", 13, 2, {0.15500512137732353, 0.014722799265496274,
		   2.9058229144393617, -55.393307878090532,
		   2.9058229144393617, -55.393307878090532,
		   133.10966174607341, -3626.8981880862489,
		   133.10966174607341, -3626.8981880862489,
		   132.35440917315739, -3604.2406108987684}},
};

const EvalueParametersByName dnaParametersByName[] = {
  {"AT77", 15, 2, {0.3361723487085515, 0.18977234047431241,
		   0.93917271154892723, -2.8454252194267946,
		   0.93917271154892723, -2.8454252194267946,
		   5.6690708632232027, -53.877468239759928,
		   5.6690708632232027, -53.877468239759928,
		   5.6546212402462572, -53.386181058543784}},

  {"ATMAP", 24, 6, {0.2308433496003553, 0.37475332645251713,
		    0.22571617180589737, -0.015715938682799857,
		    0.22571617180589737, -0.015715938682799857,
		    0.14281472963993982, 2.6440420415859895,
		    0.14281472963993982, 2.6440420415859895,
		    0.1428251000187582, 2.6434198188568869}},

  {"HOXD70", 400, 30, {0.0094463275895916698, 0.095059004086289992,
		       0.039191223120447122, -12.884975655572644,
		       0.039191223120447122, -12.884975655572644,
		       0.42791648613048661, -287.04882906895347,
		       0.42791648613048661, -287.04882906895347,
		       0.41553850600489461, -276.40376616094437}},
};

const EvalueParametersByScore dnaParametersByScore[] = {
  {1, 1, 1, 1, {0.66146689902380096, 0.026071577575191209,
		6.0609090917729898, -16.24362960699618,
		6.0609090917729898, -16.24362960699618,
		89.568020976804803, -334.27201630620215,
		89.568020976804803, -334.27201630620215,
		87.785779979559436, -327.14305231722068}},

  {1, 1, 2, 1, {0.98396641892490544, 0.17109543459934792,
		3.0348127782214962, -6.2088665291853111,
		3.0348127782214962, -6.2088665291853111,
		21.689848628115335, -94.138990367166443,
		21.689848628115335, -94.138990367166443,
		21.480894086337095, -92.885263116497001}},

  {1, 1, 7, 1, {1.0960171987681839, 0.33538787507026158,
		2.0290734315292083, -0.46514786408422282,
		2.0290734315292083, -0.46514786408422282,
		5.0543294182155085, 15.130999712620039,
		5.0543294182155085, 15.130999712620039,
		5.0543962679167036, 15.129930117400917}},

  {6, 18, 21, 9, {0.22852552944644766, 0.43337784178578154,
		  0.17722035780026291, -0.1220983494599942,
		  0.17722035780026291, -0.1220983494599942,
		  0.044639007347921783, -0.44696673676441229,
		  0.044639007347921783, -0.44696673676441229,
		  0.044598546568475124, -0.44453908999761271}},
};

const FrameshiftEvalueParameters frameshiftEvalueParameters[] = {
  {"BL62", 11, 1, 15, {0.31457181182385774, 0.077024909125411836,
		       3.5355057419386005, -39.014329998056937,
		       1.1739847579695837, -12.896364921780187,
		       160.29749789587885, -3200.2716722761552,
		       17.539096792506459, -349.06406999516196,
		       51.79266617536797, -1027.401229133852}},

  // this one is included only to make lest-test.sh faster:
  {"BL62", 11, 2, 12, {0.32704828292493776, 0.10543687626821494,
		       2.7166240938670798, -19.361170444340445,
		       0.89899140520122844, -6.29652445533966,
		       74.744763449386667, -1147.0060455603425,
		       8.3059555754127565, -127.46868078491309,
		       24.78179014970371, -379.14020451790986}},

  {"BL62", 11, 2, 15, {0.33388770870821022, 0.11532516007803961,
		       2.4678058049483518, -14.50532580281522,
		       0.82160372991753583, -4.8091552692419572,
		       55.103072639059761, -731.90592162187147,
		       6.1010321043683131, -80.763060603166991,
		       18.237750047538203, -240.59017890476582}},

  {"BL80", 11, 1, 15, {0.35400649542314511, 0.13270256108942211,
		       1.8960749679829285, -14.061923223904673,
		       0.62940827123451903, -4.6245065070665863,
		       35.18772909081801, -583.19242886423649,
		       3.8260214679558033, -62.789729751450864,
		       11.072568656496113, -178.63729131744145}},

  {"BL80", 11, 2, 15, {0.3652492341706855, 0.16850422398182774,
		       1.5316138005575799, -5.7577598061709985,
		       0.5101720323233776, -1.9097398376324572,
		       17.427364219333899, -170.0223112776693,
		       1.9259860816444827, -18.621287186644096,
		       5.7329546520583801, -54.693768145180513}},
};

static bool isEqual(const char *x, const char *y) {
  return std::strcmp(x, y) == 0;
}

static bool isHit(const EvalueParametersByName &p,
		  const char *n, int a, int b) {
  return isEqual(p.matrixName, n) && p.gapOpen == a && p.gapEpen == b;
}

static bool isHit(const EvalueParametersByScore &p,
		  int r, int q, int a, int b) {
  return p.matchScore == r && p.mismatchCost == q &&
    p.gapOpen == a && p.gapEpen == b;
}

static bool isHit(const FrameshiftEvalueParameters &p,
		  const char *n, int a, int b, int f) {
  return isEqual(p.matrixName, n) &&
    p.gapOpen == a && p.gapEpen == b && p.frameshiftCost == f;
}

static bool isProtein(const char *alphabet) {
  return isEqual(alphabet, "ACDEFGHIKLMNPQRSTVWY");
}

static bool isDna(const char *alphabet) {
  return isEqual(alphabet, "ACGT");
}

static void makeMatrix(size_t size, long *data, long **rows) {
  for (size_t i = 0; i < size; ++i)
    rows[i] = data + i * size;
}

// We have to transpose the matrix, because: for DNA-versus-protein
// alignment the input matrix is score[protein][DNA], but the ALP
// library wants score[DNA][protein].

static void copyMatrix(size_t size, const ScoreMatrixRow *in, long **out) {
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j)
      out[i][j] = in[j][i];  // transpose
}

static void copyMatrix(size_t size, const ScoreMatrixRow *in, long **out,
		       const size_t *usedLetters) {
  for (size_t i = 0; i < size; ++i)
    for (size_t j = 0; j < size; ++j)
      out[i][j] = in[usedLetters[j]][usedLetters[i]];  // transpose
}

static void makeCodonTable(long *codonTable, const GeneticCode &geneticCode) {
  unsigned char c[3];
  for (c[0] = 0; c[0] < 4; ++c[0])
    for (c[1] = 0; c[1] < 4; ++c[1])
      for (c[2] = 0; c[2] < 4; ++c[2])
	*codonTable++ = geneticCode.translation(c);
}

static size_t findUsedLetters(size_t *usedLetters, size_t alphabetSize,
			      const long *codonTable) {
  const long *end = codonTable + 64;
  size_t j = 0;
  for (size_t i = 0; i < scoreMatrixRowSize; ++i)
    if (i < alphabetSize || std::find(codonTable, end, i) < end)
      usedLetters[j++] = i;
  return j;  // returns the number of used letters
}

static void shrinkCodonTable(long *codonTable, const size_t *usedLettersBeg,
			     const size_t *usedLettersEnd) {
  for (size_t i = 0; i < 64; ++i)
    codonTable[i] = std::find(usedLettersBeg,
			      usedLettersEnd, codonTable[i]) - usedLettersBeg;
}

void LastEvaluer::init(const char *matrixName,
		       int matchScore,
		       int mismatchCost,
		       const char *alphabet,
		       const ScoreMatrixRow *scoreMatrix,
		       const double *letterFreqs1,
		       const double *letterFreqs2,
		       bool isGapped,
		       int delOpen,
		       int delEpen,
		       int insOpen,
		       int insEpen,
		       int frameshiftCost,
		       const GeneticCode &geneticCode,
		       bool isStandardGeneticCode) {
  const double lambdaTolerance = 0.01;
  const double kTolerance = 0.05;
  const double maxMegabytes = 500;
  const long randomSeed = 1;

  const double maxSeconds = 60.0;  // seems to work nicely on my computer
  size_t alphabetSize = std::strlen(alphabet);

  if (frameshiftCost > 0) {  // DNA-versus-protein alignment with frameshifts:
    if (isGapped && insOpen == delOpen && insEpen == delEpen) {
      if (isProtein(alphabet) && isStandardGeneticCode) {
	for (size_t i = 0; i < COUNTOF(frameshiftEvalueParameters); ++i) {
	  const FrameshiftEvalueParameters &p = frameshiftEvalueParameters[i];
	  if (isHit(p, matrixName, delOpen, delEpen, frameshiftCost))
	    return frameshiftEvaluer.initParameters(p.parameters);
	}
      }
    }

    long codonTable[64];
    makeCodonTable(codonTable, geneticCode);
    size_t usedLetters[scoreMatrixRowSize];
    size_t matrixSize = findUsedLetters(usedLetters, alphabetSize, codonTable);
    shrinkCodonTable(codonTable, usedLetters, usedLetters + matrixSize);

    std::vector<long> matrixData(matrixSize * matrixSize);
    std::vector<long*> matrix(matrixSize);
    makeMatrix(matrixSize, &matrixData[0], &matrix[0]);
    copyMatrix(matrixSize, scoreMatrix, &matrix[0], usedLetters);

    std::vector<double> ntFreqs(4, 0.25);
    std::vector<double> aaFreqs(matrixSize);
    copy(letterFreqs1, letterFreqs1 + alphabetSize, aaFreqs.begin());

    if (isGapped) {
      frameshiftEvaluer.initFrameshift(4, matrixSize, codonTable,
				       &matrix[0], &ntFreqs[0], &aaFreqs[0],
				       delOpen, delEpen, insOpen, insEpen,
				       frameshiftCost, true,
				       lambdaTolerance, kTolerance,
				       0, maxMegabytes, randomSeed);
    } else {
      frameshiftEvaluer.initGapless(4, matrixSize, codonTable,
				    &matrix[0], &ntFreqs[0], &aaFreqs[0]);
    }
  } else {  // ordinary alignment:
    if (isGapped && insOpen == delOpen && insEpen == delEpen) {
      if (isProtein(alphabet)) {
	for (size_t i = 0; i < COUNTOF(proteinParameters); ++i) {
	  const EvalueParametersByName &p = proteinParameters[i];
	  if (isHit(p, matrixName, delOpen, delEpen))
	    return evaluer.initParameters(p.parameters);
	}
      }
      if (isDna(alphabet)) {
	if (*matrixName) {
	  for (size_t i = 0; i < COUNTOF(dnaParametersByName); ++i) {
	    const EvalueParametersByName &p = dnaParametersByName[i];
	    if (isHit(p, matrixName, delOpen, delEpen))
	      return evaluer.initParameters(p.parameters);
	  }
	} else {
	  for (size_t i = 0; i < COUNTOF(dnaParametersByScore); ++i) {
	    const EvalueParametersByScore &p = dnaParametersByScore[i];
	    if (isHit(p, matchScore, mismatchCost, delOpen, delEpen))
	      return evaluer.initParameters(p.parameters);
	  }
	}
      }
    }

    std::vector<long> matrixData(alphabetSize * alphabetSize);
    std::vector<long*> matrix(alphabetSize);
    makeMatrix(alphabetSize, &matrixData[0], &matrix[0]);
    copyMatrix(alphabetSize, scoreMatrix, &matrix[0]);

    if (isGapped) {
      evaluer.set_gapped_computation_parameters_simplified(maxSeconds);
      evaluer.initGapped(alphabetSize, &matrix[0], letterFreqs2, letterFreqs1,
			 delOpen, delEpen, insOpen, insEpen,
			 true, lambdaTolerance, kTolerance,
			 0, maxMegabytes, randomSeed);
    } else {
      evaluer.initGapless(alphabetSize, &matrix[0],
			  letterFreqs2, letterFreqs1);
    }
  }
}

static int theMinScore(double lambda, double k, double evalue, double area) {
  // do log of evalue separately, to reduce the risk of overflow:
  double s = (std::log(k * area) - std::log(evalue)) / lambda;
  return std::ceil(std::max(1.0, s));
}

int LastEvaluer::minScore(double evalue, double area) const {
  if (evaluer.isGood()) {
    const Sls::ALP_set_of_parameters &p = evaluer.parameters();
    return theMinScore(p.lambda, p.K, evalue, area);
  } else if (frameshiftEvaluer.isGood()) {
    const Sls::FALP_set_of_parameters &p = frameshiftEvaluer.parameters();
    return theMinScore(p.lambda, p.K, evalue, area);
  } else {
    return -1;
  }
}

void LastEvaluer::writeCommented(std::ostream& out) const {
  if (evaluer.isGood()) {
    const Sls::ALP_set_of_parameters &p = evaluer.parameters();
    out << "# lambda=" << p.lambda << " K=" << p.K << "\n";
  } else if (frameshiftEvaluer.isGood()) {
    const Sls::FALP_set_of_parameters &p = frameshiftEvaluer.parameters();
    out << "# lambda=" << p.lambda << " K=" << p.K << "\n";
  }
}

void LastEvaluer::writeParameters(std::ostream &out) const {
  std::streamsize prec = out.precision(17);
  if (evaluer.isGood()) {
    const Sls::ALP_set_of_parameters &p = evaluer.parameters();
    out << p.lambda << ", " << p.K << ",\n"
	<< p.a_I << ", " << p.b_I << ",\n"
	<< p.a_J << ", " << p.b_J << ",\n"
	<< p.alpha_I << ", " << p.beta_I << ",\n"
	<< p.alpha_J << ", " << p.beta_J << ",\n"
	<< p.sigma << ", " << p.tau << "\n";
  } else if (frameshiftEvaluer.isGood()) {
    const Sls::FALP_set_of_parameters &p = frameshiftEvaluer.parameters();
    out << p.lambda << ", " << p.K << ",\n"
	<< p.a_I << ", " << p.b_I << ",\n"
	<< p.a_J << ", " << p.b_J << ",\n"
	<< p.alpha_I << ", " << p.beta_I << ",\n"
	<< p.alpha_J << ", " << p.beta_J << ",\n"
	<< p.sigma << ", " << p.tau << "\n";
  }
  out.precision(prec);
}

}
