function [ UpwindedTermX ] = matCalculateUpwindedFphysDerivative(obj, mesh, fphys, physClass, variableX )
%> @brief Function to calculate the physical variable related partial derivative in a upwind manner
%> @details
%> Function to calculate the physical variable related partial derivative in a upwind manner
%> @param[in] mesh The mesh object
%> @param[in] fphys The fphys field
%> @param[in] physClass The physical object set up
%> @param[in] termX The varialbe used to calculate the partial derivative with respect to x
%> @param[out] UpwindedTermX The partial derivative with respect to x in upwind manner

InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[vfmx, vfpx] = InnerEdge.matEvaluateSurfValue(variableX);


[ fluxX ] = obj.matEvaluateUpwindNumFlux( mesh.status, InnerEdge.FToE, ...
    fm(:,:,2), fp(:,:,2), vfmx, vfpx, ...
    InnerEdge.nx );

[ fluxMx ] = obj.matEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfmx, InnerEdge.nx );

[ fluxPx ] = obj.matEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfpx, InnerEdge.nx );

UpwindedTermX = InnerEdge.matEvaluateStrongFormEdgeRHS(fluxMx, fluxPx, fluxX);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx,...
    fm, fp, physClass.fext{1} );
[vfmx, vfpx] = BoundaryEdge.matEvaluateSurfValue(variableX);
 
[ fluxX ] = obj.matEvaluateUpwindNumFlux( mesh.status, BoundaryEdge.FToE, ...
    fm(:,:,2), fp(:,:,2), vfmx, vfpx, ...
    BoundaryEdge.nx );

[ fluxMx ] = obj.matEvaluateSurfFlux( mesh.status, BoundaryEdge.FToE, ...
    vfmx, BoundaryEdge.nx );

% the boundary edge contribution
UpwindedTermX = -UpwindedTermX - BoundaryEdge.matEvaluateStrongFormEdgeRHS(fluxMx, fluxX );

[ VolumeX ] = obj.matVolumeIntegral( mesh, cell2mat(variableX) ); 
UpwindedTermX = UpwindedTermX + VolumeX;
end