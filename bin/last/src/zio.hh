// Copyright 2017 Martin C. Frith

#ifndef ZIO_HH
#define ZIO_HH

#include "mcf_zstream.hh"

#include <iostream>
#include <iterator>  // istreambuf_iterator
#include <stdexcept>
#include <string>

namespace cbrc {

inline bool isSingleDash(const char *text) {
  return text[0] == '-' && text[1] == 0;
}

inline void openOrThrow(mcf::izstream &z, const char *fileName) {
  z.open(fileName);
  if (!z) {
    throw std::runtime_error(std::string("can't open file: ") + fileName);
  }
}

// open an input file, but if the name is "-", just return cin
inline std::istream &openIn(const char *fileName, mcf::izstream &z) {
  if (isSingleDash(fileName)) return std::cin;
  openOrThrow(z, fileName);
  return z;
}

// read a file into a string, but if the name is "-", read cin
inline std::string slurp(const char *fileName) {
  mcf::izstream z;
  std::istream &s = openIn(fileName, z);
  std::istreambuf_iterator<char> beg(s), end;
  return std::string(beg, end);
}

}

#endif
