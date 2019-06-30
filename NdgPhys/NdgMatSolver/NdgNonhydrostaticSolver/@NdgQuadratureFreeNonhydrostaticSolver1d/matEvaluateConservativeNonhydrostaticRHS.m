function RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, PhysClass)

% InnerEdge = PhysClass.meshUnion(1).InnerEdge;
% BoundaryEdge = PhysClass.meshUnion(1).BoundaryEdge;
% 
% [ fhux, ~ ] = obj.matCalculateConservativeVariableRHSMatrix( PhysClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 2);
% [ ~, fhvy ] = obj.matCalculateConservativeVariableRHSMatrix( PhysClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 3);

RHS = 2 * fphys{1}(:,:,5) + fphys{1}(:,:,1) .* obj.hux  - fphys{1}(:,:,2) .* obj.HBx;

end