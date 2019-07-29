function fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys)
mesh = physClass.meshUnion(1);
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
NonhydrostaticPressure = zeros(mesh.cell.Np, mesh.K);

Index = find(mesh.status == enumSWERegion.Wet);
NonhydroPre = reshape(NonhydroPre, mesh.cell.Np, numel(Index));
for i = 1:numel(Index)
    NonhydrostaticPressure(:,Index(i)) = NonhydroPre(:,i);
end

[ qx, qy ]  = obj.matCalculateLDGAuxialaryVariable( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydrostaticPressure,[1 2]));  % abstract in nonhydrostatic solver in 2d and 1d, different for Gauss quad and quad free verison


% fphys{1}(:,:,7) = fphys{1}(:,:,7) + NonhydrostaticPressure;
fphys{1}(:,:,7) =  NonhydrostaticPressure;
fphys{1}(:,:,6) = fphys{1}(:,:,6) + 2 * obj.dt .* NonhydrostaticPressure;
fphys{1}(:,:,2) = fphys{1}(:,:,2) - obj.dt * ( fphys{1}(:,:,1) .* qx + NonhydrostaticPressure .* obj.HBx );
fphys{1}(:,:,3) = fphys{1}(:,:,3) - obj.dt * ( fphys{1}(:,:,1) .* qy + NonhydrostaticPressure .* obj.HBy );


end