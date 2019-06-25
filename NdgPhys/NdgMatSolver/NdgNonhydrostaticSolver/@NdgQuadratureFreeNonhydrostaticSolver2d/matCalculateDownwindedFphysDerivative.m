function [ DownwindedTermX, DownwindedTermY ] = matCalculateDownwindedFphysDerivative(obj, mesh, fphys, physClass, variableX, variableY)
% Note: this function is only used for the calculation of the second order
% derivative of (H+2z)
InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[vfmx, vfpx] = InnerEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = InnerEdge.matEvaluateSurfValue(variableY);

[ fluxX, fluxY ] = mxEvaluateDownwindNumFlux( mesh.status, InnerEdge.FToE, ...
    fm(:,:,2), fm(:,:,3), fp(:,:,2), fp(:,:,3), vfmx, vfpx, vfmy, vfpy, ... 
    InnerEdge.nx, InnerEdge.ny);

[ fluxMx, fluxMy ] = mxEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfmx, vfmy, InnerEdge.nx, InnerEdge.ny);

[ fluxPx, fluxPy ] = mxEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfpx, vfpy, InnerEdge.nx, InnerEdge.ny);

% the Inner edge contribution
DownwindedTermX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxPx, fluxX);
DownwindedTermY = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMy, fluxPy, fluxY);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx,...
    BoundaryEdge.ny, fm, fp, physClass.fext{1} );
[vfmx, vfpx] = BoundaryEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = BoundaryEdge.matEvaluateSurfValue(variableY);

% vfpx = obj.matImposeNonhydroRelatedBoundaryCondition( vfmx, vfpx, enumNonhydroBoundaryCondition.ZeroGrad, obj.EidBoundaryType);
% vfpy = obj.matImposeNonhydroRelatedBoundaryCondition( vfmy, vfpy, enumNonhydroBoundaryCondition.ZeroGrad, obj.EidBoundaryType);

[ fluxX, fluxY ] = mxEvaluateDownwindNumFlux( mesh.status, BoundaryEdge.FToE, ...
    fm(:,:,2), fm(:,:,3), fp(:,:,2), fp(:,:,3), vfmx, vfpx, vfmy, vfpy, ... 
    InnerEdge.nx, InnerEdge.ny);

[ fluxMx, fluxMy ] = mxEvaluateSurfFlux( mesh.status, BoundaryEdge.FToE, ...
    vfmx, vfmy, BoundaryEdge.nx, BoundaryEdge.ny);

% the boundary edge contribution
DownwindedTermX = -DownwindedTermX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxX);
DownwindedTermY = -DownwindedTermY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMy, fluxY);

[VolumeX, VolumeY] = obj.matVolumeIntegral( mesh, cell2mat(variableX), cell2mat(variableY));
DownwindedTermX = DownwindedTermX + VolumeX;
DownwindedTermY = DownwindedTermY + VolumeY;
end