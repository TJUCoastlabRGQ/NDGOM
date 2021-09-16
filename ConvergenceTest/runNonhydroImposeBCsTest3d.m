function runNonhydroImposeBCsTest3d
% M = [ 1 0.5 0.25 0.125 ];
Mz = [3 3 3 3 ];
M = [ 2 4 8 16 ];
Order = [1 2];
Nmesh = numel(M);
Ndeg = numel(Order);
for n = 1:Ndeg
    for m = 1:Nmesh
        Solver = NonhydroImposeBCsTest3d(Order(n), Order(n), M(m), Mz(m));
        Solver.matImposeNonhydroBCsTest;
        clear Solver;
    end
end
end