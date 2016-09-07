// Copyright 2008, 2009, 2010, 2011, 2013, 2014 Martin C. Frith

// This struct holds the command line arguments for lastdb.

#ifndef LASTDB_ARGUMENTS_HH
#define LASTDB_ARGUMENTS_HH

#include <stddef.h>  // size_t
#include "SequenceFormat.hh"

#include <string>
#include <vector>

namespace cbrc{

struct LastdbArguments{
  typedef unsigned indexT;

  // set the parameters to their default values:
  LastdbArguments();

  // set parameters from a list of arguments:
  void fromArgs( int argc, char** argv, bool isOptionsOnly = false );

  // set parameters from a command line (by splitting it into arguments):
  void fromLine( const std::string& line );

  // set parameters from lines beginning with "#lastdb":
  void fromString( const std::string& s );

  void resetCumulativeOptions() { seedPatterns.clear(); verbosity = 0; }

  // options:
  bool isProtein;
  bool isKeepLowercase;
  int tantanSetting;
  bool isCaseSensitive;
  std::vector< std::string > seedPatterns;
  size_t volumeSize;  // type?
  indexT indexStep;
  std::string subsetSeedFile;
  std::string userAlphabet;
  indexT minSeedLimit;
  indexT bucketDepth;
  int childTableType;
  bool isCountsOnly;
  int verbosity;
  sequenceFormat::Enum inputFormat;

  // positional arguments:
  std::string lastdbName;
  int inputStart;  // index in argv of first input filename
};

}  // end namespace
#endif
