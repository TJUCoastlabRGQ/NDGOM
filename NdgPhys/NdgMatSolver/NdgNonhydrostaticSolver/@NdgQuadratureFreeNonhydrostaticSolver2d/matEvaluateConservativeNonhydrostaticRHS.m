function RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, PhysClass)

InnerEdge = PhysClass.meshUnion(1).InnerEdge;
BoundaryEdge = PhysClass.meshUnion(1).BoundaryEdge;

[ fhx, fhy ] = obj.matCalculateConservativeVariableRHSMatrix( PhysClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 1);
[ fhux, ~ ] = obj.matCalculateConservativeVariableRHSMatrix( PhysClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 2);
[ ~, fhvy ] = obj.matCalculateConservativeVariableRHSMatrix( PhysClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 3);

RHS = -2 * fphys{1}(:,:,6) + 2 * fphys{1}(:,:,2) .* obj.bx  + 2 * fphys{1}(:,:,3) .* obj.by + ...
  fhx.* fphys{1}(:,:,2) + fhy.* fphys{1}(:,:,3)  - fhux.* fphys{1}(:,:,1) - fhvy.* fphys{1}(:,:,1);
end