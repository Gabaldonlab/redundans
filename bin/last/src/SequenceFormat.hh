// Copyright 2011 Martin C. Frith

#ifndef SEQUENCE_FORMAT_HH
#define SEQUENCE_FORMAT_HH

#include <cctype>
#include <istream>
#include <stddef.h>
#include <string>

namespace cbrc {

namespace sequenceFormat {
enum Enum { fastx, fastqSanger, fastqSolexa, fastqIllumina, prb, pssm, fasta,
	    fastxKeep };
}

inline std::istream &operator>>(std::istream &s, sequenceFormat::Enum &f) {
  std::string w;
  s >> w;
  if (!s) return s;
  for (size_t i = 0; i < w.size(); ++i) {
    w[i] = std::tolower(w[i]);
  }
  /**/ if (w == "0" || w == "fastx") f = sequenceFormat::fastx;
  else if (w == "keep") f = sequenceFormat::fastxKeep;
  else if (w == "1" || w == "sanger") f = sequenceFormat::fastqSanger;
  else if (w == "2" || w == "solexa") f = sequenceFormat::fastqSolexa;
  else if (w == "3" || w == "illumina") f = sequenceFormat::fastqIllumina;
  else if (w == "4" || w == "prb") f = sequenceFormat::prb;
  else if (w == "5" || w == "pssm") f = sequenceFormat::pssm;
  else s.setstate(std::ios::failbit);
  return s;
}

inline bool isUseFastq(sequenceFormat::Enum f) {
  return
      f == sequenceFormat::fastqSanger ||
      f == sequenceFormat::fastqSolexa ||
      f == sequenceFormat::fastqIllumina;
}

inline bool isUseQuality(sequenceFormat::Enum f) {
  return isUseFastq(f) || f == sequenceFormat::prb;
}

inline int qualityOffset(sequenceFormat::Enum f) {
  return (f == sequenceFormat::fastqSanger) ? 33 : 64;
}  // The result is meaningless for non-quality formats.

inline bool isPhred(sequenceFormat::Enum f) {
  return
      f == sequenceFormat::fastqSanger || f == sequenceFormat::fastqIllumina;
}

}

#endif
