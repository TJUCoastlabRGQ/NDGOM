function [ DownwindedTermX, DownwindedTermY ] = matCalculateDownwindedFphysDerivative(obj, mesh, fphys, physClass, variableX, variableY)
% Note: this function is only used for the calculation of the second order
% derivative of (H+2z)
InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[vfmx, vfpx] = InnerEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = InnerEdge.matEvaluateSurfValue(variableY);

% getCentralFluxTerm
fluxX = InnerEdge.nx .* ( vfmx + vfpx )./2;
fluxY = InnerEdge.ny .* ( vfmy + vfpy )./2;

% getUpwindFluxTerm
temphum = fm(:,:,2); temphvm = fm(:,:,3);
temphup = fp(:,:,2); temphvp = fp(:,:,3);
Index = ( temphum .* InnerEdge.nx + temphvm .* InnerEdge.ny > 0 &...
    - temphup .* InnerEdge.nx - temphvp .* InnerEdge.ny <= 0 );
fluxX(Index) = vfpx(Index) .* InnerEdge.nx(Index);
fluxY(Index) = vfpy(Index) .* InnerEdge.ny(Index);

Index = ( temphum .* InnerEdge.nx + temphvm .* InnerEdge.ny <= 0 & ...
    - temphup .* InnerEdge.nx - temphvp .* InnerEdge.ny > 0 );
fluxX( Index ) =  vfmx(Index) .* InnerEdge.nx(Index);
fluxY(Index) = vfmy(Index) .* InnerEdge.ny(Index);
% wet-dry interface considered
[vfmx, vfpx, fluxX] = obj.matGetWetDryFaceVariableAndFlux( vfmx, vfpx, fluxX );
[vfmy, vfpy, fluxY] = obj.matGetWetDryFaceVariableAndFlux( vfmy, vfpy, fluxY );

% getLocalAndAdjacentFluxTerm
fluxMx = InnerEdge.nx .* vfmx; fluxMy = InnerEdge.ny .* vfmy;
fluxPx = InnerEdge.nx .* vfpx; fluxPy = InnerEdge.ny .* vfpy;

% the Inner edge contribution
DownwindedTermX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxPx, fluxX);
DownwindedTermY = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMy, fluxPy, fluxY);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx,...
    BoundaryEdge.ny, fm, fp, physClass.fext{1} );
[vfmx, vfpx] = BoundaryEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = BoundaryEdge.matEvaluateSurfValue(variableY);

vfpx = obj.matImposeNonhydroRelatedBoundaryCondition( vfmx, vfpx, enumNonhydroBoundaryCondition.ZeroGrad, obj.EidBoundaryType);
vfpy = obj.matImposeNonhydroRelatedBoundaryCondition( vfmy, vfpy, enumNonhydroBoundaryCondition.ZeroGrad, obj.EidBoundaryType);

% getCentralFluxTerm
fluxX = BoundaryEdge.nx .* ( vfmx + vfpx )./2; fluxY = BoundaryEdge.ny .* ( vfmy + vfpy )./2;
temphum = fm(:,:,2); temphvm = fm(:,:,3);
temphup = fp(:,:,2); temphvp = fp(:,:,3);

Index = ( temphum .* BoundaryEdge.nx + temphvm .* BoundaryEdge.ny > 0 & ...
    - temphup .* BoundaryEdge.nx - temphvp .* BoundaryEdge.ny <= 0 );
fluxX(Index) = vfpx(Index) .* BoundaryEdge.nx(Index);
fluxY(Index) = vfpy(Index) .* BoundaryEdge.ny(Index);

Index = ( temphum .* BoundaryEdge.nx + temphvm .* BoundaryEdge.ny <= 0 & ...
    - temphup .* BoundaryEdge.nx - temphvp .* BoundaryEdge.ny > 0 );
fluxX( Index ) =  vfmx(Index) .* BoundaryEdge.nx(Index);
fluxY( Index ) =  vfmy(Index) .* BoundaryEdge.ny(Index);

% getLocalFluxTerm
fluxMx = BoundaryEdge.nx .* vfmx; fluxMy = BoundaryEdge.ny .* vfmy;

% the boundary edge contribution
DownwindedTermX = -DownwindedTermX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxX);
DownwindedTermY = -DownwindedTermY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMy, fluxY);

[VolumeX, VolumeY] = obj.matVolumeIntegral( mesh, cell2mat(variableX), cell2mat(variableY));
DownwindedTermX = DownwindedTermX + VolumeX;
DownwindedTermY = DownwindedTermY + VolumeY;
end