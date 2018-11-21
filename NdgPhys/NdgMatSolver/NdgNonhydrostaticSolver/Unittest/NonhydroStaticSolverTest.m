Norder = [1 2 3];
% for i = 1:numel(Norder)
%     Solver = NdgFlatBottomNonhydrostaticTest(Norder(i), NdgCellType.Quad);
% %     Solver.testmatEvaluateNonhydrostaticSurfaceValue;
% %     Solver.testmatEvaluateVelocityLikeTerms;
% %     Solver.testmatEvaluatePenaltyTerm;
% %     Solver.testmatEvaluateNonconservativeFlux;
% %     Solver.testmatEvaluateDeltaSurfaceFlux;
%     Solver.testmatEvaluateDivergenceDirectionX;
%     Solver.testmatEvaluateDivergenceDirectionY;
%     Solver.testmatEvaluateDeltaSurfaceFluxX;
%     Solver.testmatEvaluateDeltaSurfaceFluxY;
%   
% end

for i = 1:numel(Norder)
    Solver = NdgFlatBottomWetDryTestWithClampedBoundary(Norder(i), NdgCellType.Quad);
    Solver.testAssembleWetDryInterface;
    Solver.testEidBoundaryType;
    Solver.testReverseEidBoundaryType;
    Solver.testGetFaceOuterValue;
    Solver.testFluxTerm;
    clear Solver;
end

for i = 1:numel(Norder)
    Solver = NdgSlopeBottomWetDryTestWithWallBoundary(Norder(i), NdgCellType.Quad);
    Solver.testAssembleWetDryInterface;
    Solver.testTopographyNonhydrostaticFaceOuterValue;
    clear Solver;
end

for i = 1:numel(Norder)
    Solver = NdgFlatBottomWetDryTestWithWallBoundary(Norder(i), NdgCellType.Quad);
    Solver.testAssembleWetDryInterface;
    Solver.testEidBoundaryType;
    Solver.testReverseEidBoundaryType;
    Solver.testUpdatedExactEidBoundaryType;
    Solver.testGetFaceOuterValue;
    clear Solver;  
end