// Author: Martin C. Frith 2022
// SPDX-License-Identifier: GPL-3.0-or-later

#ifndef LAST_SPLIT_OPTIONS_HH
#define LAST_SPLIT_OPTIONS_HH

#include <stddef.h>

#include <string>
#include <vector>

struct LastSplitOptions {
  int format;
  bool isTopSeqQuery;
  std::string genome;
  int direction;
  double cis;
  double trans;
  double mean;
  double sdev;
  double mismap;
  int score;
  bool no_split;
  size_t bytes;
  bool verbose;
  bool isSplicedAlignment;
  std::vector<std::string> inputFileNames;

  static const char helpf[];
  static const char helpr[];
  static const char helpg[];
  static const char helpd[];
  static const char helpc[];
  static const char helpt[];
  static const char helpM[];
  static const char helpS[];
  static const char helpm[];
  static const char helps[];
  static const char helpn[];
  static const char helpb[];

  LastSplitOptions();

  void setUnspecifiedValues(int lastalMinScore, double scale);

  void print() const;

  static char parseOutputFormat(const char *text);
};

#endif
