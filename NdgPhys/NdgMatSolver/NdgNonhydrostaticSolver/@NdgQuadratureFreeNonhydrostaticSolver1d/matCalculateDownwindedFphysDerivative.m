function [ DownwindedTermX ] = matCalculateDownwindedFphysDerivative(obj, mesh, fphys, physClass, variableX)
% Note: this function is only used for the calculation of the second order
% derivative of (H+2z)
InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[vfmx, vfpx] = InnerEdge.matEvaluateSurfValue(variableX);

[ fluxX ] = mxEvaluateDownwindNumFlux( mesh.status, InnerEdge.FToE, ...
    fm(:,:,2), fp(:,:,2), vfmx, vfpx, InnerEdge.nx);

[ fluxMx ] = mxEvaluateSurfFlux( mesh.status, InnerEdge.FToE, vfmx, InnerEdge.nx );

[ fluxPx ] = mxEvaluateSurfFlux( mesh.status, InnerEdge.FToE, ...
    vfpx, InnerEdge.nx );

% fluxMx = InnerEdge.nx .* vfmx;
% fluxPx = InnerEdge.nx .* vfpx;

% the Inner edge contribution
DownwindedTermX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxPx, fluxX);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx,...
    fm, fp, physClass.fext{1} );

[vfmx, vfpx] = BoundaryEdge.matEvaluateSurfValue(variableX);

vfpx = obj.matImposeNonhydroRelatedBoundaryCondition( vfmx, vfpx, enumNonhydroBoundaryCondition.ZeroGrad, obj.EidBoundaryType);

[ fluxX ] = mxEvaluateDownwindNumFlux( mesh.status, BoundaryEdge.FToE, ...
    fm(:,:,2), fp(:,:,2), vfmx, vfpx, BoundaryEdge.nx );

[ fluxMx ] = mxEvaluateSurfFlux( mesh.status, BoundaryEdge.FToE, vfmx, BoundaryEdge.nx );
% fluxMx = BoundaryEdge.nx .* vfmx;
% the boundary edge contribution
DownwindedTermX = -DownwindedTermX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxX);

[ VolumeX ] = obj.matVolumeIntegral( mesh, cell2mat(variableX));
DownwindedTermX = DownwindedTermX + VolumeX;
end