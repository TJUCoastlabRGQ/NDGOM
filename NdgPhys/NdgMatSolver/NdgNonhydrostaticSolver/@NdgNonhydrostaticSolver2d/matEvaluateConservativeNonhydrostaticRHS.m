function RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass)

InnerEdge = physClass.meshUnion(1).InnerEdge;
BoundaryEdge = physClass.meshUnion(1).BoundaryEdge;
fhx = obj.matCalculateConservativeVariableRelatedMatrixX(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 1);
fhy = obj.matCalculateConservativeVariableRelatedMatrixY(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 1);
fhux = obj.matCalculateConservativeVariableRelatedMatrixX(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 2);
fhvy = obj.matCalculateConservativeVariableRelatedMatrixY(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 3);
RHS = -2 * ((fphys{1}(:,:,6) - fphys{1}(:,:,2) .* obj.bx  - fphys{1}(:,:,3) .* obj.by)) + ...
     ((fhx.* fphys{1}(:,:,2))) + ((fhy.* fphys{1}(:,:,3)))  - ((fhux.* fphys{1}(:,:,1))) - ...
     ((fhvy.* fphys{1}(:,:,1)));
 
end