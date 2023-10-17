// Copyright 2016 Martin C. Frith

#ifndef THREAD_UTIL_HH
#define THREAD_UTIL_HH

#include <iostream>
#include <stdexcept>

#ifdef HAS_CXX_THREADS
#include <thread>
#endif

namespace cbrc {

inline unsigned decideNumberOfThreads(unsigned request,
				      const char *programName,
				      bool isVerbose) {
#ifdef HAS_CXX_THREADS
  if (request) return request;
  unsigned x = std::thread::hardware_concurrency();
  if (x) {
    if (isVerbose) std::cerr << programName << ": threads=" << x << '\n';
    return x;
  }
  std::cerr
    << programName
    << ": can't determine how many threads to use: falling back to 1 thread\n";
#else
  if (request != 1)
    throw
      std::runtime_error("I was installed here with multi-threading disabled");
#endif
  return 1;
}

}

#endif
