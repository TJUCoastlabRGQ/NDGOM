function [ UpwindedTermX, UpwindedTermY ] = matCalculateUpwindedFphysDerivative(obj, mesh, fphys, physClass, variableX, variableY)
%> @brief Function to calculate the physical variable related partial derivative in a upwind manner
%> @details
%> Function to calculate the physical variable related partial derivative in a upwind manner
%> @param[in] mesh The mesh object
%> @param[in] fphys The fphys field
%> @param[in] physClass The physical object set up
%> @param[in] termX The varialbe used to calculate the partial derivative with respect to x
%> @param[in] termY The varialbe used to calculate the partial derivative with respect to y
%> @param[out] UpwindedTermX The partial derivative with respect to x in upwind manner
%> @param[out] UpwindedTermY The partial derivative with respect to y in upwind manner

InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[vfmx, vfpx] = InnerEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = InnerEdge.matEvaluateSurfValue(variableY);


[ fluxX, fluxY ] = obj.matEvaluateUpwindNumFlux( mesh.status, InnerEdge.FToE, ...
    fm(:,:,2), fm(:,:,3), fp(:,:,2), fp(:,:,3), vfmx, vfpx, vfmy, vfpy, ...
    InnerEdge.nx, InnerEdge.ny);

% [ fluxMx, fluxMy ] = mxEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
%     vfmx, vfmy, InnerEdge.nx, InnerEdge.ny);
[ fluxMx, fluxMy ] = obj.matEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfmx, vfmy, InnerEdge.nx, InnerEdge.ny);

% [ fluxPx, fluxPy ] = mxEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
%     vfpx, vfpy, InnerEdge.nx, InnerEdge.ny);
[ fluxPx, fluxPy ] = obj.matEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfpx, vfpy, InnerEdge.nx, InnerEdge.ny);

UpwindedTermX = InnerEdge.matEvaluateStrongFormEdgeRHS(fluxMx, fluxPx, fluxX);
UpwindedTermY = InnerEdge.matEvaluateStrongFormEdgeRHS(fluxMy, fluxPy, fluxY);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx,...
    BoundaryEdge.ny, fm, fp, physClass.fext{1} );
[vfmx, vfpx] = BoundaryEdge.matEvaluateSurfValue(variableX);
[vfmy, vfpy] = BoundaryEdge.matEvaluateSurfValue(variableY);

% [ fluxX, fluxY ] = mxEvaluateUpwindNumFlux( mesh.status, BoundaryEdge.FToE, ...
%     fm(:,:,2), fm(:,:,3), fp(:,:,2), fp(:,:,3), vfmx, vfpx, vfmy, vfpy, ...
%     BoundaryEdge.nx, BoundaryEdge.ny);
[ fluxX, fluxY ] = obj.matEvaluateUpwindNumFlux( mesh.status, BoundaryEdge.FToE, ...
    fm(:,:,2), fm(:,:,3), fp(:,:,2), fp(:,:,3), vfmx, vfpx, vfmy, vfpy, ...
    BoundaryEdge.nx, BoundaryEdge.ny);

% [ fluxMx, fluxMy ] = mxEvaluateSurfFlux( mesh.status, BoundaryEdge.FToE, ...
%     vfmx, vfmy, BoundaryEdge.nx, BoundaryEdge.ny);
[ fluxMx, fluxMy ] = obj.matEvaluateSurfFlux( mesh.status, BoundaryEdge.FToE, ...
    vfmx, vfmy, BoundaryEdge.nx, BoundaryEdge.ny);

% the boundary edge contribution
UpwindedTermX = -UpwindedTermX - BoundaryEdge.matEvaluateStrongFormEdgeRHS(fluxMx, fluxX );
UpwindedTermY = -UpwindedTermY - BoundaryEdge.matEvaluateStrongFormEdgeRHS(fluxMy, fluxY );

[VolumeX, VolumeY] = obj.matVolumeIntegral( mesh, cell2mat(variableX), cell2mat(variableY));   %abstract in nonhydrostatic solver 2d and 1d, different for quad free and quad version
UpwindedTermX = UpwindedTermX + VolumeX;
UpwindedTermY = UpwindedTermY + VolumeY;
end