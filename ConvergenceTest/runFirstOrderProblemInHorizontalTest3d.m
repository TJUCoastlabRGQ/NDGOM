function runFirstOrderProblemInHorizontalTest3d
% M = [ 1 0.5 0.25 0.125 ];
% Mz = [10 20 40 80 ];
M = [ 2 4 8 16 ];
Order = [1 2 3];
Nmesh = numel(M);
Ndeg = numel(Order);
for n = 1:Ndeg
    for m = 1:Nmesh
        Solver = FirstOrderProblemInHorizontalTest3d(Order(n), Order(n), M(m));
        Solver.EllipticProblemSolve;
        clear Solver;
    end
end
end