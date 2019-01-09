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

% NqHx = matCalculateUpwindedNonhydroRelatedMatrixX(obj, physClass, BoundaryEdge, InnerEdge,  num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);
% NqHy = matCalculateUpwindedNonhydroRelatedMatrixY(obj, physClass, BoundaryEdge, InnerEdge,  num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);

% [NqHx, NqHy] = obj.matEvaluateLocalDerivativeTerm(mesh, NonhydroVolumeflux);

fphys{1}(:,:,6) = fphys{1}(:,:,6) + obj.dt/physClass.rho .* NonhydrostaticPressure;
fphys{1}(:,:,2) = fphys{1}(:,:,2) - obj.dt/physClass.rho .* (NqHx +NonhydrostaticPressure.*obj.bx);
fphys{1}(:,:,3) = fphys{1}(:,:,3) - obj.dt/physClass.rho .* (NqHy +NonhydrostaticPressure.*obj.by);

% UpdataWaterDepth(physClass, mesh, fphys, tempFluxhu, tempFluxhv);

end

function termX = matCalculateUpwindedNonhydroRelatedMatrixX(obj, physClass, BoundaryEdge, InnerEdge, Variable, ftype)
%> @brief Function to calculate the characteristic matrix
%> @details Function to calculate the characteristic matrix in the x direction
%> @param[in] BoundaryEdge the boundary edge object
%> @param[in] InnerEdge the inner edge object
%> @param[in] Variable variable used to calculate the characteristic matrix
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[out] termX the calculated characteristic matrix in x direction
%> @param[out] termY the calculated characteristic matrix in y direction
%< Inner value and outer value of the Inner edges
mesh = physClass.meshUnion(1);
[fm, fp] = InnerEdge.matEvaluateSurfValue( Variable );       
[fm, fp] = obj.matGetFaceValue(fm, fp, ftype);
%< Inner edge contribution
fluxMX = InnerEdge.nx.*fm;
fluxPX = InnerEdge.nx.*fp; 
% fluxSx = InnerEdge.nx .* (1 + sign(InnerEdge.nx .* fm))./2 .* fm + InnerEdge.nx .*  (1 + sign( -InnerEdge.nx .* fp))./2 .* fp; % upwind
fluxSx = InnerEdge.nx .* (1 + sign(InnerEdge.nx .* fm))./2 .* fp + InnerEdge.nx .*  (1 + sign( -InnerEdge.nx .* fp))./2 .* fm; % downwind
termX = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMX, fluxPX, fluxSx );

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( Variable );        
fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, ftype, obj.EidBoundaryType);
%% test first
%< Boundary edge contribution
fluxMX = BoundaryEdge.nx.*fm;
% fluxSX = BoundaryEdge.ny .* (1 + sign(BoundaryEdge.ny .* fm))./2 .* fm+ BoundaryEdge.ny .*  (1 + sign( -BoundaryEdge.ny .* fp))./2 .* fp;
fluxSX = BoundaryEdge.ny .* (1 + sign(BoundaryEdge.ny .* fm))./2 .* fp+ BoundaryEdge.ny .*  (1 + sign( -BoundaryEdge.ny .* fp))./2 .* fm;
termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

termX = termX + mesh.rx .* (mesh.cell.Dr * cell2mat(Variable))...
    + mesh.sx .* (mesh.cell.Ds * cell2mat(Variable));

end

function termY = matCalculateUpwindedNonhydroRelatedMatrixY(obj, physClass, BoundaryEdge, InnerEdge, Variable, ftype)
%> @brief Function to calculate the characteristic matrix
%> @details Function to calculate the characteristic matrix in the x direction
%> @param[in] BoundaryEdge the boundary edge object
%> @param[in] InnerEdge the inner edge object
%> @param[in] Variable variable used to calculate the characteristic matrix
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[out] termX the calculated characteristic matrix in x direction
%> @param[out] termY the calculated characteristic matrix in y direction
%< Inner value and outer value of the Inner edges
mesh = physClass.meshUnion(1);
[fm, fp] = InnerEdge.matEvaluateSurfValue( Variable );       
[fm, fp] = obj.matGetFaceValue(fm, fp, ftype);
%< Inner edge contribution
fluxMY = InnerEdge.ny.*fm;
fluxPY = InnerEdge.ny.*fp; 
% fluxSy = InnerEdge.ny .* (1 + sign(InnerEdge.ny .* fm))./2 .* fm + InnerEdge.ny .*  (1 + sign( -InnerEdge.ny .* fp))./2 .* fp;
fluxSy = InnerEdge.ny .* (1 + sign(InnerEdge.ny .* fm))./2 .* fp + InnerEdge.ny .*  (1 + sign( -InnerEdge.ny .* fp))./2 .* fm;
termY = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMY, fluxPY, fluxSy );

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( Variable );        
fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, ftype, obj.EidBoundaryType);
%% test first
%< Boundary edge contribution
fluxMY = BoundaryEdge.ny.*fm;
% fluxSy = BoundaryEdge.ny .* (1 + sign(BoundaryEdge.ny .* fm))./2 .* fm +  BoundaryEdge.ny .*  (1 + sign( -BoundaryEdge.ny .* fp))./2 .* fp;  % upwind
fluxSy = BoundaryEdge.ny .* (1 + sign(BoundaryEdge.ny .* fm))./2 .* fp + BoundaryEdge.ny .*  (1 + sign( -BoundaryEdge.ny .* fp))./2 .* fm;  % downwind
termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSy);

termY = termY + mesh.ry .* (mesh.cell.Dr * cell2mat(Variable))...
    + mesh.sy .* (mesh.cell.Ds * cell2mat(Variable));

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