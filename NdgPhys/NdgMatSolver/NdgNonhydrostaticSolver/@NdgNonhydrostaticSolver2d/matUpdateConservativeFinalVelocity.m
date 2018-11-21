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

NonhydroVolumeflux = 1/2 * NonhydrostaticPressure .* fphys{1}(:,:,1);

NqHx = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);
NqHy = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);

fphys{1}(:,:,6) = fphys{1}(:,:,6) + obj.dt/physClass.rho  .* NonhydrostaticPressure;
fphys{1}(:,:,2) = fphys{1}(:,:,2) + (- obj.dt/physClass.rho) *(NqHx +NonhydrostaticPressure.*obj.bx);
fphys{1}(:,:,3) = fphys{1}(:,:,3) + (- obj.dt/physClass.rho) *(NqHy +NonhydrostaticPressure.*obj.by);

% UpdataWaterDepth(physClass, mesh, fphys, tempFluxhu, tempFluxhv);

end

function  UpdataWaterDepth(physClass, mesh, fphys, tempFluxhu, tempFluxhv)

dt = physClass.dt;
resQ = zeros(size(fphys{1}(:,:,1)));
deltaFluxhu = fphys{1}(:,:,2) - tempFluxhu;  deltaFluxhv = fphys{1}(:,:,3) - tempFluxhv;
% deltaFluxx = (deltaFluxhu(mesh.eidM) + deltaFluxhu(mesh.eidP))/2;
% deltaFluxy = (deltaFluxhv(mesh.eidM) + deltaFluxhv(mesh.eidP))/2;
deltaFluxx = ((1 + sign(deltaFluxhu(mesh.eidM))).*deltaFluxhu(mesh.eidM) +  (1 - sign(deltaFluxhu(mesh.eidP))).*deltaFluxhu(mesh.eidP))./2;
deltaFluxy = ((1 + sign(deltaFluxhv(mesh.eidM))).*deltaFluxhv(mesh.eidM) +  (1 - sign(deltaFluxhv(mesh.eidP))).*deltaFluxhv(mesh.eidP))./2;
index = (mesh.eidtype ~= 0);
deltaFluxx(index) = 0;deltaFluxy(index) = 0;

fluxq = mesh.nx .* ( deltaFluxhu(mesh.eidM) - deltaFluxx) + mesh.ny .* ( deltaFluxhv(mesh.eidM) - deltaFluxy);

[rk4a, rk4b, ~] = GetRKParamter();

for intRK = 1:5
    
    RHS = GetRhs(mesh, deltaFluxhu, deltaFluxhv, fluxq);
    
    resQ = rk4a(intRK)*resQ + dt*RHS;
    fphys{1}(:,:,1) ...
        = fphys{1}(:,:,1) + rk4b(intRK)*resQ;
    
end
end

function rhs = GetRhs(mesh, deltaFluxhu, deltaFluxhv, fluxq)
rhs =  - (mesh.rx .* (mesh.cell.Dr * deltaFluxhu) + mesh.sx .* (mesh.cell.Ds * deltaFluxhu) + ...
   mesh.ry .* (mesh.cell.Dr * deltaFluxhv) + mesh.sy .* (mesh.cell.Ds * deltaFluxhv) - mesh.cell.LIFT*(mesh.Js .* fluxq)./mesh.J);
end


function [rk4a, rk4b, rk4c] = GetRKParamter()
rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0];
end