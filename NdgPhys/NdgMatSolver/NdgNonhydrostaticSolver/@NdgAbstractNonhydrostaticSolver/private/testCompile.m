% ipath = ['-I','D:\����\Sparse�������\laspack\laspack' ];
% mex('-v', '-largeArrayDims', ipath, 'mxEvaluateSparseNonhydrostaticPressure.c', '\����\Sparse�������\laspack\laspack\matrix.c')

% mex('-v', '-largeArrayDims', 'mxEvaluateSparseNonhydrostaticPressure.c', 'laspack\matrix.c')

%  mex -v -largeArrayDims mxEvaluateSparseNonhydrostaticPressure.c matrix.c laspack/errhandl.c

%doable %  mex -g -v -largeArrayDims mxEvaluateSparseNonhydrostaticPressure.c laspack/matrix.c errhandl.c    

%doable %mex -g -largeArrayDims mxEvaluateSparseNonhydrostaticPressure.c laspack/matrix.c laspack/errhandl.c  %quick

mex -g -largeArrayDims mxEvaluateSparseNonhydrostaticPressure.c laspack/*.c      %time slow