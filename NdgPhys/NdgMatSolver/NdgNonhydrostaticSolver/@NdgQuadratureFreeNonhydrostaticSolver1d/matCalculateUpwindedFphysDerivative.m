function [ UpwindedTermX ] = matCalculateUpwindedFphysDerivative(obj, mesh, fphys, physClass, variableX )

InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[vfmx, vfpx] = InnerEdge.matEvaluateSurfValue(variableX);

[ fluxX ] = mxEvaluateUpwindNumFlux( mesh.status, InnerEdge.FToE, ...
    fm(:,:,2), fp(:,:,2), vfmx, vfpx, InnerEdge.nx );

% getLocalAndAdjacentFluxTerm
fluxMx = InnerEdge.nx .* vfmx;
fluxPx = InnerEdge.nx .* vfpx;

% the Inner edge contribution
UpwindedTermX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxPx, fluxX);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx,...
     fm, fp, physClass.fext{1} );
[vfmx, vfpx] = BoundaryEdge.matEvaluateSurfValue(variableX);

[ fluxX ] = mxEvaluateUpwindNumFlux( mesh.status, BoundaryEdge.FToE, ...
    fm(:,:,2), fp(:,:,2), vfmx, vfpx, BoundaryEdge.nx );

[ fluxMx ] = mxEvaluateSurfFlux( mesh.status, BoundaryEdge.FToE, vfmx, BoundaryEdge.nx );

% the boundary edge contribution
UpwindedTermX = -UpwindedTermX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxX);

[ VolumeX ] = obj.matVolumeIntegral( mesh, cell2mat(variableX));
UpwindedTermX = UpwindedTermX + VolumeX;
end