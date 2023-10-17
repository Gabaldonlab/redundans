// Author: Martin C. Frith 2017
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef GETOPT_UTIL_HH
#define GETOPT_UTIL_HH

#include <getopt.h>

inline void resetGetopt() {  // XXX fragile voodoo
#ifdef __GLIBC__
  optind = 0;
#else
  optind = 1;
  //optreset = 1;  // XXX ???
#endif
}

#endif
