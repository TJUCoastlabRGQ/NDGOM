function runAssembleFinalStiffMatrixTest3d
%RUNASSEMBLEFINALSTIFFMATRIXTEST3D 此处显示有关此函数的摘要
%   此处显示详细说明
Mz = [3 3 3 3 ];
M = [ 2 4 8 16 ];
Order = [1 2];
Nmesh = numel(M);
Ndeg = numel(Order);
for n = 1:Ndeg
    for m = 1:Nmesh
        Solver = EllipticProblemAssembleFinalStiffMatrix3d(Order(n), Order(n), M(m), Mz(m));
        Solver.matStiffMatrixCompare;
        clear Solver;
    end
end
end

