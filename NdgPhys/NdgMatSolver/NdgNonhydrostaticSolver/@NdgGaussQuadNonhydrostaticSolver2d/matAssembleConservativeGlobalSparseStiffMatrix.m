function StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix( obj, PhysClass, UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY, UpdatedNP, fphys)

mesh = PhysClass.meshUnion(1);
[ Hx, Hy ] = obj.matCalculateConservativeVariableRHSMatrix( PhysClass, mesh, fphys, enumNonhydroBoundaryCondition.Zero, 1);
Qx = UpdatedPNPX; Qy = UpdatedPNPY;
Q2x = UpdatedSPNPX; Q2y = UpdatedSPNPY;
H = fphys{1}(:,:,1);
for i = 1:numel(fphys{1}(:,:,1))
    Qx(:,i) = H(i) .* Qx(:,i);
    Qy(:,i) = H(i) .* Qy(:,i);
    Q2x(:,i) = H(i) .* Q2x(:,i);
    Q2y(:,i) = H(i) .* Q2y(:,i);    
end

StiffMatrix = UpdatedPNPX + UpdatedPNPY + UpdatedSPNPX + UpdatedSPNPY + UpdatedNP;

tempH = obj.Vq{1} * H;

for i = 1:numel(fphys{1}(:,:,1))
StiffMatrix(:,i) =  ( 2 * obj.dt * UpdatedNP(:,i) + obj.dt * ( Hx(:) .* Qx(:,i)) + obj.dt * ( Hy(:) .* Qy(:,i)) - ...
    obj.dt * ( tempH(:) .* Q2x(:,i) ) - obj.dt * ( tempH(:) .* Q2y(:,i) )) ./ 1000;
end
end