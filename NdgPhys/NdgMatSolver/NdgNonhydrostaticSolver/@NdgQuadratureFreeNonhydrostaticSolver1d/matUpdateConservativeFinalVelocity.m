function fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys)
mesh = physClass.meshUnion(1);
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
NonhydrostaticPressure = zeros(mesh.cell.Np, mesh.K);
NonhydroPre = reshape(NonhydroPre, mesh.cell.Np, obj.WetNum);
for i = 1:numel(obj.WetCellIndex)
    NonhydrostaticPressure(:,obj.WetCellIndex(i)) = NonhydroPre(:,i);
end

[ qx ]  = obj.matCalculateLDGAuxialaryVariable( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydrostaticPressure,[1 2]));


% fphys{1}(:,:,6) = fphys{1}(:,:,6) + NonhydrostaticPressure;
fphys{1}(:,:,6) =  NonhydrostaticPressure;
fphys{1}(:,:,5) = fphys{1}(:,:,6) + 2 * obj.dt .* NonhydrostaticPressure;
fphys{1}(:,:,2) = fphys{1}(:,:,2) - obj.dt * ( fphys{1}(:,:,1) .* qx + NonhydrostaticPressure .* obj.HBx );
end