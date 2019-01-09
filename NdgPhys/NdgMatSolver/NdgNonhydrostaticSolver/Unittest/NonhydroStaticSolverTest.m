Norder = [1 2 3];
for i = 1:numel(Norder)
    %> This part is  mainly used to test the global matrix assemble
    Solver = NdgFlatBottomNonhydrostaticTest(Norder(i), enumStdCell.Quad);
    Solver.testGlobalMatrixAssemble;
    Solver.testEidBoundaryType;
    clear Solver;
end

for i = 1:numel(Norder)
    %> This part is  mainly used to test the non-hydrostatic related boundary value 
    Solver = NdgFlatBottomNonhydrostaticBoundaryValueTest(Norder(i), enumStdCell.Quad);
    Solver.testEidBoundaryType;
    Solver.testmatImposeNonhydroRelatedBoundaryCondition;
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

% What to do next: change the matGetFaceValue, get and test the fm and fp value, for this part, 
% we plan to take both the fm and fp at the wet-dry interface to be zero at the zero condition boundary