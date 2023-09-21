% minFunc
fprintf('Compiling minFunc files...\n');
%mex minFunc/mcholC.c
%mex minFunc/lbfgsC.c 

mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims minFunc/mcholC.c
mex CFLAGS="\$CFLAGS -std=c99" -largeArrayDims minFunc/lbfgsC.c