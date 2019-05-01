function RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, PhysClass)

% InnerEdge = PhysClass.meshUnion(1).InnerEdge;
% BoundaryEdge = PhysClass.meshUnion(1).BoundaryEdge;
mesh = PhysClass.meshUnion(1);
[ fhx, fhy ] = obj.matCalculateConservativeVariableRHSMatrix( PhysClass, mesh, fphys, enumNonhydroBoundaryCondition.Zero, 1);
[ fhux, ~ ] =  obj.matCalculateConservativeVariableRHSMatrix( PhysClass, mesh, fphys, enumNonhydroBoundaryCondition.Zero, 2);
[ ~, fhvy ] =  obj.matCalculateConservativeVariableRHSMatrix( PhysClass, mesh, fphys, enumNonhydroBoundaryCondition.Zero, 3);

HW = obj.Vq{1} * fphys{1}(:,:,6);
HU = obj.Vq{1} * fphys{1}(:,:,2);
HV = obj.Vq{1} * fphys{1}(:,:,3);
H  = obj.Vq{1} * fphys{1}(:,:,1);
% RHS = -2 * fphys{1}(:,:,6) + 2 * fphys{1}(:,:,2) .* obj.bx  + 2 * fphys{1}(:,:,3) .* obj.by + ...
%   fhx.* fphys{1}(:,:,2) + fhy.* fphys{1}(:,:,3)  - fhux.* fphys{1}(:,:,1) - fhvy.* fphys{1}(:,:,1);
RHS = -2 * HW + fhx.* HU + fhy.* HV  - fhux.* H - fhvy.* H;
end