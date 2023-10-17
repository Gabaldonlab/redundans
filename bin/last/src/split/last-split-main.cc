// Author: Martin C. Frith 2013
// SPDX-License-Identifier: GPL-3.0-or-later

// Parse the command line and run last-split.

#include "last-split.hh"

#include "stringify.hh"

#include <getopt.h>

#include <cstdlib>  // EXIT_SUCCESS, EXIT_FAILURE
#include <iostream>

static void run(int argc, char* argv[]) {
  LastSplitOptions opts;

  std::string version = "last-split "
#include "version.hh"
    "\n";

  std::string help = "\
Usage: " + std::string(argv[0]) + " [options] LAST-alignments.maf\n\
\n\
Read alignments of query sequences to a genome, and estimate the genomic\n\
source of each part of each query, allowing different parts of one query to\n\
come from different parts of the genome.\n\
\n\
Options:\n\
 -h, --help         show this help message and exit\n\
 -f, --format=FMT   " + LastSplitOptions::helpf + "\n\
 -r, --reverse      " + LastSplitOptions::helpr + "\n\
 -g, --genome=NAME  " + LastSplitOptions::helpg + "\n\
 -d, --direction=D  " + LastSplitOptions::helpd + "\n\
 -c, --cis=PROB     " + LastSplitOptions::helpc + "\n\
 -t, --trans=PROB   " + LastSplitOptions::helpt + "\n\
 -M, --mean=MEAN    " + LastSplitOptions::helpM + "\n\
 -S, --sdev=SDEV    " + LastSplitOptions::helpS + "\n\
 -m, --mismap=PROB  " + LastSplitOptions::helpm + "\n\
 -s, --score=INT    " + LastSplitOptions::helps + "\n\
 -n, --no-split     " + LastSplitOptions::helpn + "\n\
 -b, --bytes=B      " + LastSplitOptions::helpb + "\n\
 -v, --verbose      be verbose\n\
 -V, --version      show version information and exit\n\
";

  const char sOpts[] = "hf:rg:d:c:t:M:S:m:s:nb:vV";

  static struct option lOpts[] = {
    { "help",     no_argument,       0, 'h' },
    { "format",   required_argument, 0, 'f' },
    { "reverse",  no_argument,       0, 'r' },
    { "genome",   required_argument, 0, 'g' },
    { "direction",required_argument, 0, 'd' },
    { "cis",      required_argument, 0, 'c' },
    { "trans",    required_argument, 0, 't' },
    { "mean",     required_argument, 0, 'M' },
    { "sdev",     required_argument, 0, 'S' },
    { "mismap",   required_argument, 0, 'm' },
    { "score",    required_argument, 0, 's' },
    { "no-split", no_argument,       0, 'n' },
    { "bytes",    required_argument, 0, 'b' },
    { "verbose",  no_argument,       0, 'v' },
    { "version",  no_argument,       0, 'V' },
    { 0, 0, 0, 0}
  };

  int c;
  while ((c = getopt_long(argc, argv, sOpts, lOpts, &c)) != -1) {
    switch (c) {
    case 'h':
      std::cout << help;
      return;
    case 'f':
      opts.format = LastSplitOptions::parseOutputFormat(optarg);
      break;
    case 'r':
      opts.isTopSeqQuery = true;
      break;
    case 'g':
      opts.isSplicedAlignment = true;
      opts.genome = optarg;
      break;
    case 'd':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.direction, optarg);
      break;
    case 'c':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.cis, optarg);
      break;
    case 't':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.trans, optarg);
      break;
    case 'M':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.mean, optarg);
      break;
    case 'S':
      opts.isSplicedAlignment = true;
      cbrc::unstringify(opts.sdev, optarg);
      break;
    case 'm':
      cbrc::unstringify(opts.mismap, optarg);
      break;
    case 's':
      cbrc::unstringify(opts.score, optarg);
      break;
    case 'n':
      opts.no_split = true;
      break;
    case 'b':
      cbrc::unstringifySize(opts.bytes, optarg);
      break;
    case 'v':
      opts.verbose = true;
      break;
    case 'V':
      std::cout << version;
      return;
    case '?':
      throw std::runtime_error("");
    }
  }

  opts.inputFileNames.assign(argv + optind, argv + argc);

  if (opts.inputFileNames.empty()) opts.inputFileNames.push_back("-");

  std::ios_base::sync_with_stdio(false);  // makes std::cin much faster!!!

  lastSplit(opts);
}

int main(int argc, char* argv[]) {
  try {
    run(argc, argv);
    if (!flush(std::cout)) throw std::runtime_error("write error");
    return EXIT_SUCCESS;
  } catch (const std::bad_alloc& e) {  // bad_alloc::what() may be unfriendly
    std::cerr << argv[0] << ": out of memory\n";
    return EXIT_FAILURE;
  } catch (const std::exception& e) {
    const char *s = e.what();
    if (*s) std::cerr << argv[0] << ": " << s << '\n';
    return EXIT_FAILURE;
  }
}
