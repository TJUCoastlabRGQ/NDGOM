function [ UpwindedTermX, UpwindedTermY ] = matCalculateUpwindedFphysDerivative(obj, mesh, fphys, physClass, variableX, variableY)

InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[vfmx, vfpx] = InnerEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = InnerEdge.matEvaluateSurfValue(variableY);

[ fluxX, fluxY ] = mxEvaluateUpwindNumFlux( mesh.status, InnerEdge.FToE, ...
    fm(:,:,2), fm(:,:,3), fp(:,:,2), fp(:,:,3), vfmx, vfpx, vfmy, vfpy, ... 
    InnerEdge.nx, InnerEdge.ny);

[ fluxMx, fluxMy ] = mxEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfmx, vfmy, InnerEdge.nx, InnerEdge.ny);

[ fluxPx, fluxPy ] = mxEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfpx, vfpy, InnerEdge.nx, InnerEdge.ny);

UpwindedTermX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxPx, fluxX);
UpwindedTermY = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMy, fluxPy, fluxY);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx,...
    BoundaryEdge.ny, fm, fp, physClass.fext{1} );
[vfmx, vfpx] = BoundaryEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = BoundaryEdge.matEvaluateSurfValue(variableY);

[ fluxX, fluxY ] = mxEvaluateUpwindNumFlux( mesh.status, BoundaryEdge.FToE, ...
    fm(:,:,2), fm(:,:,3), fp(:,:,2), fp(:,:,3), vfmx, vfpx, vfmy, vfpy, ... 
    BoundaryEdge.nx, BoundaryEdge.ny);

[ fluxMx, fluxMy ] = mxEvaluateSurfFlux( mesh.status, BoundaryEdge.FToE, ...
    vfmx, vfmy, BoundaryEdge.nx, BoundaryEdge.ny);

% the boundary edge contribution
UpwindedTermX = -UpwindedTermX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxX );
UpwindedTermY = -UpwindedTermY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMy, fluxY );

[VolumeX, VolumeY] = obj.matVolumeIntegral( mesh, cell2mat(variableX), cell2mat(variableY));
UpwindedTermX = UpwindedTermX + VolumeX;
UpwindedTermY = UpwindedTermY + VolumeY;
end