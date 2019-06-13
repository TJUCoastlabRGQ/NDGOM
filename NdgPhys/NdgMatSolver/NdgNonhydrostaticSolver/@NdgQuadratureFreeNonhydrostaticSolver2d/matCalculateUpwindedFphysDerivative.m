function [ UpwindedTermX, UpwindedTermY ] = matCalculateUpwindedFphysDerivative(obj, mesh, fphys, physClass, variableX, variableY)

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
fluxX(Index) = vfmx(Index) .* InnerEdge.nx(Index);
fluxY(Index) = vfmy(Index) .* InnerEdge.ny(Index);

Index = ( temphum .* InnerEdge.nx + temphvm .* InnerEdge.ny <= 0 & ...
    - temphup .* InnerEdge.nx - temphvp .* InnerEdge.ny > 0 );
fluxX( Index ) =  vfpx(Index) .* InnerEdge.nx(Index);
fluxY(Index) = vfpy(Index) .* InnerEdge.ny(Index);
% wet-dry interface considered
[vfmx, vfpx, fluxX] = obj.matGetWetDryFaceVariableAndFlux( vfmx, vfpx, fluxX );
[vfmy, vfpy, fluxY] = obj.matGetWetDryFaceVariableAndFlux( vfmy, vfpy, fluxY );

% getLocalAndAdjacentFluxTerm
fluxMx = InnerEdge.nx .* vfmx; fluxMy = InnerEdge.ny .* vfmy;
fluxPx = InnerEdge.nx .* vfpx; fluxPy = InnerEdge.ny .* vfpy;

% the Inner edge contribution
UpwindedTermX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxPx, fluxX);
UpwindedTermY = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMy, fluxPy, fluxY);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx,...
    BoundaryEdge.ny, fm, fp, physClass.fext{1} );
[vfmx, vfpx] = BoundaryEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = BoundaryEdge.matEvaluateSurfValue(variableY);

% getCentralFluxTerm
fluxX = BoundaryEdge.nx .* ( vfmx + vfpx )./2; fluxY = BoundaryEdge.ny .* ( vfmy + vfpy )./2;
temphum = fm(:,:,2); temphvm = fm(:,:,3);
temphup = fp(:,:,2); temphvp = fp(:,:,3);

% getUpwindFluxTerm
Index = ( temphum .* BoundaryEdge.nx + temphvm .* BoundaryEdge.ny > 0 & ...
    - temphup .* BoundaryEdge.nx - temphvp .* BoundaryEdge.ny <= 0 );
fluxX(Index) = vfmx(Index) .* BoundaryEdge.nx(Index);
fluxY(Index) = vfmy(Index) .* BoundaryEdge.ny(Index);

Index = ( temphum .* BoundaryEdge.nx + temphvm .* BoundaryEdge.ny <= 0 & ...
    - temphup .* BoundaryEdge.nx - temphvp .* BoundaryEdge.ny > 0 );
fluxX( Index ) =  vfpx(Index) .* BoundaryEdge.nx(Index);
fluxY( Index ) =  vfpy(Index) .* BoundaryEdge.ny(Index);

% getLocalFluxTerm
fluxMx = BoundaryEdge.nx .* vfmx; fluxMy = BoundaryEdge.ny .* vfmy;

% the boundary edge contribution
UpwindedTermX = -UpwindedTermX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxX );
UpwindedTermY = -UpwindedTermY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMy, fluxY );

[VolumeX, VolumeY] = obj.matVolumeIntegral( mesh, cell2mat(variableX), cell2mat(variableY));
UpwindedTermX = UpwindedTermX + VolumeX;
UpwindedTermY = UpwindedTermY + VolumeY;
end