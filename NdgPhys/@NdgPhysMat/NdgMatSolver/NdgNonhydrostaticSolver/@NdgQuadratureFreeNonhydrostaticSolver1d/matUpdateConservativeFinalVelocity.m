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

[ qx ]  = obj.matCalculateIPDGAuxialaryVariable( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydrostaticPressure,[1 2])); 

% fphys{1}(:,:,6) = fphys{1}(:,:,6) + NonhydrostaticPressure;
fphys{1}(:,:,6) =  NonhydrostaticPressure;
fphys{1}(:,:,5) = fphys{1}(:,:,5) + 2 * obj.dt .* NonhydrostaticPressure;
fphys{1}(:,:,2) = fphys{1}(:,:,2) - obj.dt * ( fphys{1}(:,:,1) .* qx + NonhydrostaticPressure .* obj.HBx );

end