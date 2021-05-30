function runEllipticProblemMatrixAssembleTest3dNew
%RUNELLIPTICPROBLEMTEST3DNEW 此处显示有关此函数的摘要
%   此处显示详细说明
clear, clc
M = [2.5 1.25 0.625 0.625/2 0.625/4];
Mz = [4 8 16 32 64];
Order = [1 2];
for n = 1:Ndeg
    for m = 1:Nmesh
        Solver = EllipticProblemMatrixAssembleTest3dNew(Order(n), Order(n), M(m), Mz(m));
        Solver.EllipticProblemSolve;
    end
end
end