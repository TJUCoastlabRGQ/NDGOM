function fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys)
mesh = physClass.meshUnion(1);
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
% NonhydrostaticPressure = reshape(NonhydrostaticPressure, mesh.cell.Np, mesh.K);
NonhydrostaticPressure = zeros(mesh.cell.Np, mesh.K);
NonhydroPre = reshape(NonhydroPre, mesh.cell.Np, obj.WetNum);
for i = 1:numel(obj.WetCellIndex)
    NonhydrostaticPressure(:,obj.WetCellIndex(i)) = NonhydroPre(:,i);
end

% NonhydrostaticPressure(:,1:17) = - NonhydrostaticPressure(:,1:17);
% index = ( abs(NonhydrostaticPressure) < 1 );
% NonhydrostaticPressure(index) = 0;

NonhydroVolumeflux = 1/2 * NonhydrostaticPressure .* fphys{1}(:,:,1);

% NqHx = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);
% NqHy = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);

[ NqHx , NqHy ] = obj.matCalculateCharacteristicMatrix( mesh,  BoundaryEdge, InnerEdge, num2cell(NonhydroVolumeflux,[1 2]), num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);

% NqHx = matCalculateUpwindedNonhydroRelatedMatrixX(obj, physClass, BoundaryEdge, InnerEdge,  NonhydroVolumeflux, fphys, enumNonhydroBoundaryCondition.Zero);
% NqHy = matCalculateUpwindedNonhydroRelatedMatrixY(obj, physClass, BoundaryEdge, InnerEdge,  NonhydroVolumeflux, fphys, enumNonhydroBoundaryCondition.Zero);

% data = fphys{1}(:,1,:);

fphys{1}(:,:,6) = fphys{1}(:,:,6) + obj.dt/physClass.rho .* NonhydrostaticPressure;
fphys{1}(:,:,2) = fphys{1}(:,:,2) - obj.dt/physClass.rho .* (NqHx +NonhydrostaticPressure.*obj.bx);
fphys{1}(:,:,3) = fphys{1}(:,:,3) - obj.dt/physClass.rho .* (NqHy +NonhydrostaticPressure.*obj.by);

% fphys{1}(:,1,:) = data( :,1,:);


end