function  [qx, qy, q2x, q2y, qbx, qby, fqbx, fqby] = matAssembleCharacteristicMatrix(obj, mesh, index)
%<@brief function used to calculate the term used when assemble the global stiff matrix
%<@detail In this version, the zero boundary condition for nonhydrostatic pressure is imposed at
%< the general clamped boundary, and the zero grad boundary condition at the wall boundary.
%< Also the zero gradient boundary condition for the auxialary variable is imposed at the
%< general clamped boundary, and the zero boundary condition at the wall boundary
%< @param[in] mesh The mesh object
%< @param[in] index The number of the studied point
%< @param[out] qx first derivative of nonhydrostatic pressure with respect to x
%< @param[out] qy first derivative of nonhydrostatic pressure with respect to y
%< @param[out] q2x second derivative of nonhydrostatic pressure with respect to x
%< @param[out] q2y second derivative of nonhydrostatic pressure with respect to y
%< @param[out] qbx product of nonhydrostatic pressure with topography gradient bx in x direction
%< @param[out] qby product of nonhydrostatic pressure with topography gradient by in y direction
%< @param[out] fqbx product of nonhydrostatic pressure with topography gradient bx in x direction, this term is further included in the flux term
%< @param[out] fqby product of nonhydrostatic pressure with topography gradient by in y direction, this term is further included in the flux term
K = mesh.K; Np = mesh.cell.Np;
Nonhydro = zeros( Np, K );
gmat = zeros( Np, K );
gmat(index) = 1.0/2.0;
Nonhydro(index) = 1;
qbx = zeros(K*Np,1);qby = zeros(K*Np,1);

InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;

ele = ceil(index/Np);
[UpWindedFlag, DownWindedFlag] = AssembleWindedFlagInformation(InnerEdge, ele);

[tempqx, tempqy] = obj.matCalculateFluxUpwindedTerm( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.Zero);
% tempqx = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
% tempqy = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
qx = tempqx(:);
qy = tempqy(:);
% [nqx, nqy] = obj.matEvaluateLocalDerivativeTerm( mesh, gmat );
% qx = nqx(:);
% qy = nqy(:);

[tempq2x, ~] = obj.matCalculateFluxDownwindedTerm( mesh, BoundaryEdge, InnerEdge, num2cell(tempqx,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.ZeroGrad);
[~, tempq2y] = obj.matCalculateFluxDownwindedTerm( mesh, BoundaryEdge, InnerEdge, num2cell(tempqy,[1 2]),UpWindedFlag, DownWindedFlag, enumNonhydroBoundaryCondition.ZeroGrad);
% tempq2x = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(tempqx,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
% tempq2y = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(tempqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
% tempq2x = matCalculatePenaltyCharacteristicMatrixX(obj, mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), num2cell(tempqx,[1 2]), obj.EidBoundaryType);
% tempq2y = matCalculatePenaltyCharacteristicMatrixY(obj, mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), num2cell(tempqy,[1 2]), obj.EidBoundaryType);
q2x = tempq2x(:); q2y = tempq2y(:);

% tempfqbx = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(Nonhydro.*obj.bx,[1 2]), enumNonhydroBoundaryCondition.Zero);
% tempfqby = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(Nonhydro.*obj.by,[1 2]), enumNonhydroBoundaryCondition.Zero);
tempfqbx = zeros(size(mesh.x));
tempfqby = zeros(size(mesh.x));
fqbx = tempfqbx(:);
fqby = tempfqby(:);
qbx(index) = obj.bx(index);qby(index) = obj.by(index);
end

function [UpWindedFlag, DownWindedFlag] = AssembleWindedFlagInformation(InnerEdge, ele)
DownWindedFlag = zeros(size(InnerEdge.nx));
UpWindedFlag = ones(size(InnerEdge.nx));
[Row, Col] = find(InnerEdge.FToE == ele);
for i = 1:numel(Row)
    if Row(i) == 2
        UpWindedFlag(:,Col(i)) = 0;
        DownWindedFlag(:, Col(i)) = 1;
    end
end
end

function tempq2x = matCalculatePenaltyCharacteristicMatrixX(obj, mesh, BoundaryEdge, InnerEdge, gmat, tempqx, EidBoundaryType)
tau = 1.0;
%% Inner Edge part
[um, up] = InnerEdge.matEvaluateSurfValue(gmat);
[um, up] = obj.matGetFaceValue(um, up, enumNonhydroBoundaryCondition.Zero);
JumpU = InnerEdge.nx .* (um - up) + InnerEdge.ny .* (um - up);   % Calculate the jump term finished

[qm, qp] = InnerEdge.matEvaluateSurfValue(tempqx);
[qm, qp] = obj.matGetFaceValue(qm, qp, enumNonhydroBoundaryCondition.ZeroGrad);
fluxq = (qm + qp)./2.*InnerEdge.nx - tau * JumpU.*InnerEdge.nx; %Calculate the numerical flux finished

fluxM = InnerEdge.nx .* qm;   
fluxP = InnerEdge.nx .* qp;

tempq2x = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxM, fluxP, fluxq); %Inner Edge Integration finished

%% Boundary Edge part
[um, up] = BoundaryEdge.matEvaluateSurfValue(gmat);
up = obj.matImposeNonhydroRelatedBoundaryCondition(um, up, enumNonhydroBoundaryCondition.Zero, EidBoundaryType);
JumpU = BoundaryEdge.nx .* (um - up) + BoundaryEdge.ny .* (um - up);   % Calculate the jump term finished


[qm, qp] = BoundaryEdge.matEvaluateSurfValue(tempqx);
qp = obj.matImposeNonhydroRelatedBoundaryCondition(qm, qp, enumNonhydroBoundaryCondition.ZeroGrad, EidBoundaryType);
fluxq = (qm + qp)./2.*BoundaryEdge.nx - tau * JumpU.*BoundaryEdge.nx; %Calculate the numerical flux finished

fluxM = BoundaryEdge.nx .* qm; 

tempq2x = -tempq2x - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxM, fluxq);

tempq2x = tempq2x + mesh.rx .* (mesh.cell.Dr * cell2mat(tempqx))...
    + mesh.sx .* (mesh.cell.Ds * cell2mat(tempqx));
end

function tempq2y = matCalculatePenaltyCharacteristicMatrixY(obj, mesh, BoundaryEdge, InnerEdge, gmat, tempqy, EidBoundaryType )
tau = 10.0;
%% Inner Edge part
[um, up] = InnerEdge.matEvaluateSurfValue(gmat);
[um, up] = obj.matGetFaceValue(um, up, enumNonhydroBoundaryCondition.Zero);
JumpU = InnerEdge.nx .* (um - up) + InnerEdge.ny .* (um - up);   % Calculate the jump term finished

[qm, qp] = InnerEdge.matEvaluateSurfValue(tempqy);
[qm, qp] = obj.matGetFaceValue(qm, qp, enumNonhydroBoundaryCondition.ZeroGrad);
fluxq = (qm + qp)./2.*InnerEdge.ny - tau * JumpU.*InnerEdge.ny; %Calculate the numerical flux finished

fluxM = InnerEdge.ny .* qm;   
fluxP = InnerEdge.ny .* qp;

tempq2y = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxM, fluxP, fluxq); %Inner Edge Integration finished

%% Boundary Edge part
[um, up] = BoundaryEdge.matEvaluateSurfValue(gmat);
up = obj.matImposeNonhydroRelatedBoundaryCondition(um, up, enumNonhydroBoundaryCondition.Zero, EidBoundaryType);
JumpU = BoundaryEdge.nx .* (um - up) + BoundaryEdge.ny .* (um - up);   % Calculate the jump term finished


[qm, qp] = BoundaryEdge.matEvaluateSurfValue(tempqy);
qp = obj.matImposeNonhydroRelatedBoundaryCondition(qm, qp, enumNonhydroBoundaryCondition.ZeroGrad, EidBoundaryType);
fluxq = (qm + qp)./2.*BoundaryEdge.ny - tau * JumpU.*BoundaryEdge.ny; %Calculate the numerical flux finished

fluxM = BoundaryEdge.ny .* qm; 

tempq2y = -tempq2y - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxM, fluxq);

tempq2y = tempq2y + mesh.ry .* (mesh.cell.Dr * cell2mat(tempqy))...
    + mesh.sy .* (mesh.cell.Ds * cell2mat(tempqy));
end


