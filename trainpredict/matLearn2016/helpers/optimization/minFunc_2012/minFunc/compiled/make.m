% mac 2023b


mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims ../mex/lbfgsAddC.c
mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims ../mex/lbfgsC.c
mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims ../mex/lbfgsProdC.c
mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims ../mex/mcholC.c
		