Norder = [1 2 3];
for i = 1:numel(Norder)
    Solver = NdgFlatBottomNonhydrostaticTest(Norder(i), enumStdCell.Quad);
    Solver.testGlobalMatrixAssemble;
    Solver.testEidBoundaryType;
    clear Solver;
end

for i = 1:numel(Norder)
    Solver = NdgFlatBottomWetDryTestWithClampedBoundary(Norder(i), enumStdCell.Quad);
    Solver.testAssembleWetDryInterface;
    Solver.testEidBoundaryType;
    Solver.testmatGetFaceValue;
    clear Solver;
end

for i = 1:numel(Norder)
    Solver = NdgFlatBottomWetDryTestWithWallBoundary(Norder(i), enumStdCell.Quad);
    Solver.testAssembleWetDryInterface;
    Solver.testEidBoundaryType;
    clear Solver;  
end