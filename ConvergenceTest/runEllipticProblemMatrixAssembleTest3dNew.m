function runEllipticProblemMatrixAssembleTest3dNew
%RUNELLIPTICPROBLEMTEST3DNEW �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
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