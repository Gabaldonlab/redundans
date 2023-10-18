// Copyright 2017 Martin C. Frith

// mcf::izstream is similar to std::ifstream.  The difference is, if
// you give it a gzip-compressed file, it will decompress what it
// reads.

#ifndef MCF_ZSTREAM_HH
#define MCF_ZSTREAM_HH

#include <zlib.h>

#include <cstdio>  // BUFSIZ
#include <istream>
#include <stdexcept>
#include <streambuf>

namespace mcf {

class zbuf : public std::streambuf {
public:
  zbuf() : input(0) {}

  ~zbuf() { close(); }

  bool is_open() const { return input; }

  zbuf *open(const char *fileName) {
    if (is_open()) return 0;
    input = gzopen(fileName, "rb");
    if (!is_open()) return 0;
    return this;
  }

  zbuf *close() {
    if (!is_open()) return 0;
    int e = gzclose(input);
    input = 0;
    return (e == Z_OK || e == Z_BUF_ERROR) ? this : 0;
  }

protected:
  int underflow() {
    if (gptr() == egptr()) {
      int size = gzread(input, buffer, BUFSIZ);
      if (size < 0) throw std::runtime_error("gzread error");
      setg(buffer, buffer, buffer + size);
    }
    return (gptr() == egptr()) ?
      traits_type::eof() : traits_type::to_int_type(*gptr());
  }

private:
  gzFile input;
  char buffer[BUFSIZ];
};

class izstream : public std::istream {
public:
  izstream() : std::istream(&buf) {}

  izstream(const char *fileName) : std::istream(&buf) {
    open(fileName);
  }

  bool is_open() const { return buf.is_open(); }

  void open(const char *fileName) {
    // do something special if fileName is "-"?
    if (!buf.open(fileName)) setstate(failbit);
    else clear();
  }

  void close() {
    if (!buf.close()) setstate(failbit);
  }

private:
  zbuf buf;
};

}

#endif
