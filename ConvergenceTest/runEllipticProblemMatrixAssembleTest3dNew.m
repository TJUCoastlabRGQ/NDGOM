function runEllipticProblemMatrixAssembleTest3dNew
%RUNELLIPTICPROBLEMTEST3DNEW 此处显示有关此函数的摘要
%   此处显示详细说明
clear, clc
M = [10 5 2.5 1.25];
Mz = [2 4 8 16];
Order = [1 2];

for n = 1:numel(Order)
    for m = 1:numel(M)
        Solver = EllipticProblemMatrixAssembleTest3dNew(Order(n), Order(n), M(m), Mz(m));
        Solver.EllipticProblemSolve;
    end
end
end