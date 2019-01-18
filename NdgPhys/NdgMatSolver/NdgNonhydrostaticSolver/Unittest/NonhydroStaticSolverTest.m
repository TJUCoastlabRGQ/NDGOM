% For all the test follows, enumNonhydroBoundaryCondition.Zero
% corresponds to the eidtype of boundary edge flaged to be -1, for those
% boundary edge point, the non-hydrostatic pressure there is set to be
% zero. Likewise, enumNonhydroBoundaryCondition.ZeroGrad corresponds to the
% same edge point with the non-hydrostatic pressure gradient there set to be zero


Norder = [1 2 3];
% for i = 1:numel(Norder)
%     %> This part is  mainly used to test the non-hydrostatic related boundary value 
%     Solver = NdgFlatBottomNonhydrostaticBoundaryValueTest(Norder(i), enumStdCell.Quad);
%     display('=========Test for boundary type begin==========');    
%     Solver.testEidBoundaryType;
%     display('=========Test for boundary type finished==========');   
%     Solver.testmatImposeNonhydroRelatedBoundaryCondition;
%     clear Solver;
% end
% 
% for i = 1:numel(Norder)
%     %> This part is mainly used to test matrix used to assemble the global
%     %> stiff matrix and the conservative variable related partial derivatives
%     Solver = NdgFlatBottomNonhydrostaticMatrixTest(Norder(i), enumStdCell.Quad);
%     display('=========Test for characteristic matrix in x direction begin==========');
%     Solver.testmatCalculateCharacteristicMatrixX;
%     display('=========Test for characteristic matrix in x direction finished==========');
%     display('=========Test for characteristic matrix in y direction begin==========');
%     Solver.testmatCalculateCharacteristicMatrixY;
%     display('=========Test for characteristic matrix in y direction finished==========');
%     display('=========Test for conservative variable related matrix in x direction begin==========');
%     Solver.testmatCalculateConservativeVariableRelatedMatrixX;
%     display('=========Test for conservative variable related matrix in x direction finished==========');
%     display('=========Test for conservative variable related matrix in y direction begin==========');
%     Solver.testmatCalculateConservativeVariableRelatedMatrixY;
%     display('=========Test for conservative variable related matrix in y direction finished==========');
%     display('=========Global stiff matrix test begin==========');
%     Solver.testGlobalMatrix;
%     display('=========Global stiff matrix test finished==========');    
%     clear Solver;
% end
% 
% for i = 1:numel(Norder)
%     %> This part is the supplementation of the last test. Compared with the
%     %> last one, the initial field setting is different
%     Solver = NdgFlatBottomNonhydrostaticMatrixSecondTest(Norder(i), enumStdCell.Quad);
%     display('=========Test for characteristic matrix in x direction begin==========');
%     Solver.testmatCalculateCharacteristicMatrixX;
%     display('=========Test for characteristic matrix in x direction finished==========');
%     display('=========Test for characteristic matrix in y direction begin==========');    
%     Solver.testmatCalculateCharacteristicMatrixY;
%     display('=========Test for characteristic matrix in y direction finished==========');
%     display('=========Test for conservative variable related matrix in x direction begin==========');    
%     Solver.testmatCalculateConservativeVariableRelatedMatrixX;
%     display('=========Test for conservative variable related matrix in x direction finished==========');
%     display('=========Test for conservative variable related matrix in y direction begin==========');    
%     Solver.testmatCalculateConservativeVariableRelatedMatrixY;
%     display('=========Test for conservative variable related matrix in y direction finished==========');
%     display('=========Global stiff matrix test begin==========');
%     Solver.testGlobalMatrix;
%     display('=========Global stiff matrix test finished==========');    
%     clear Solver;
% end
% 
% for i = 1:numel(Norder)
%     %> This part is  mainly used to test the global matrix assemble
%     Solver = NdgFlatBottomNonhydrostaticTest(Norder(i), enumStdCell.Quad);
%     display('=========Test for the global stiff matrix assemble begin==========');    
%     Solver.testGlobalMatrixAssemble;
%     display('=========Test for the global stiff matrix assemble finished==========');    
%     Solver.testEidBoundaryType;
%     clear Solver;
% end


% %Actually, this part has been finished in the program with central flux for
% %global water depth derivative and mixed upwinded and central flux term for the flux term
% for i = 1:numel(Norder)
%     %> This part is  mainly used to test the calculation of the upwinded advective term in the RHS
%     Solver = NdgFlatBottomNonhydrostaticUpwindedTermTest(Norder(i), enumStdCell.Quad);
%     display('=========Test for the calculation of the upwinded advective term in the RHS begin==========');    
%     Solver.testmatCalculateConservativeVariableRelatedUpwindedMatrixXAndY;
%     display('=========Test for the calculation of the upwinded advective term in the RHS finished==========');    
%     clear Solver;    
% 
% end

for i = 1:numel(Norder)
    %> This part is  mainly used to test the calculation of the upwinded advective term in the RHS
    Solver = NdgFlatBottomNonhydrostaticLDGFluxTermTest(Norder(i), enumStdCell.Quad);
    display('=========Test for the calculation of the LDG flux term for global matrix assemble begin==========');    
    Solver.testInnerEdgeScalerJumpTerm;
    Solver.testInnerEdgeVecrtorJumpTerm;
    
    Solver.testBoundaryEdgeScalerJumpTerm;
    Solver.testBoundaryEdgeVecrtorJumpTerm;
    display('=========Test for the calculation of the LDG flux term for global matrix assemble finished==========');    
    clear Solver;    
end

% for i = 1:numel(Norder)
%     Solver = NdgFlatBottomWetDryTestWithClampedBoundary(Norder(i), enumStdCell.Quad);
%     Solver.testAssembleWetDryInterface;
%     Solver.testEidBoundaryType;
%     Solver.testmatGetFaceValue;
%     clear Solver;
% end
% 
% for i = 1:numel(Norder)
%     Solver = NdgFlatBottomWetDryTestWithWallBoundary(Norder(i), enumStdCell.Quad);
%     Solver.testAssembleWetDryInterface;
%     Solver.testEidBoundaryType;
%     clear Solver;  
% end

% What to do next: change the matGetFaceValue, get and test the fm and fp value, for this part, 
% we plan to take both the fm and fp at the wet-dry interface to be zero at the zero condition boundary