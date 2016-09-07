const struct {
  const char *name;
  const char *text;
} subsetSeeds[] = {
{"BISF", "\
#lastdb -w2\n\
#lastal -pBISF\n\
1  CT A G\n\
0  ACGT\n\
1111110101100\n\
"},

{"BISR", "\
#lastdb -w2\n\
#lastal -pBISR\n\
1  AG C T\n\
0  ACGT\n\
1111110101100\n\
"},

{"MAM4", "\
1  A C G T\n\
0  ACGT\n\
T  AG CT\n\
11100TT01T00T10TTTT\n\
TTTT110TT0T001T0T1T1\n\
11TT010T01TT0001T\n\
11TT10T1T101TT\n\
"},

{"MAM8", "\
1  A C G T\n\
0  ACGT\n\
T  AG CT\n\
1101T1T0T1T00TT1TT\n\
1TTTTT010TT0TT01011TTT\n\
1TTTT10010T011T0TTTT1\n\
111T011T0T01T100\n\
1T10T100TT01000TT01TT11\n\
111T101TT000T0T10T00T1T\n\
111100T011TTT00T0TT01T\n\
1T1T10T1101101\n\
"},

{"MURPHY10", "\
#lastdb -p\n\
1  ILMV FWY A C G H P KR ST DENQ\n\
1\n\
"},

{"NEAR", "\
#lastal -r6 -q18 -a21 -b9\n\
1  A C G T\n\
0  ACGT\n\
1111110\n\
"},

{"YASS", "\
1  A C G T\n\
0  ACGT\n\
T  AG CT\n\
1T1001100101\n\
"},

};
