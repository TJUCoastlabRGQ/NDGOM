Solver = WindDrivenFlow(1,1,100,10);
Solver.dt = 5;
a = [1 2 3 4 5 6 7 8 9];
b = [9 8 7 6 5 4 3 2 1];
test(struct(Solver.meshUnion.InnerEdge), Solver.dt, int8(Solver.meshUnion.type), a, b);
disp(Solver.dt);