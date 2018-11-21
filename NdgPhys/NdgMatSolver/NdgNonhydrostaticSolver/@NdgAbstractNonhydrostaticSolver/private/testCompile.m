% ipath = ['-I','D:\资料\Sparse矩阵求解\laspack\laspack' ];
% mex('-v', '-largeArrayDims', ipath, 'mxEvaluateSparseNonhydrostaticPressure.c', '\资料\Sparse矩阵求解\laspack\laspack\matrix.c')

% mex('-v', '-largeArrayDims', 'mxEvaluateSparseNonhydrostaticPressure.c', 'laspack\matrix.c')

%  mex -v -largeArrayDims mxEvaluateSparseNonhydrostaticPressure.c matrix.c laspack/errhandl.c

%doable %  mex -g -v -largeArrayDims mxEvaluateSparseNonhydrostaticPressure.c laspack/matrix.c errhandl.c    

%doable %mex -g -largeArrayDims mxEvaluateSparseNonhydrostaticPressure.c laspack/matrix.c laspack/errhandl.c  %quick

mex -g -largeArrayDims mxEvaluateSparseNonhydrostaticPressure.c laspack/*.c      %time slow