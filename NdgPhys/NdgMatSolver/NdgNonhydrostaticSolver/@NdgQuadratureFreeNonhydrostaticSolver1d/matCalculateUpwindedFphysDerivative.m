function [ UpwindedTermX ] = matCalculateUpwindedFphysDerivative(obj, mesh, fphys, physClass, variableX )

InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[vfmx, vfpx] = InnerEdge.matEvaluateSurfValue(variableX);

% getCentralFluxTerm
fluxX = InnerEdge.nx .* ( vfmx + vfpx )./2;

% getUpwindFluxTerm
temphum = fm(:,:,2); 
temphup = fp(:,:,2); 
Index = ( temphum .* InnerEdge.nx  > 0 & - temphup .* InnerEdge.nx  <= 0 );
fluxX(Index) = vfmx(Index) .* InnerEdge.nx(Index);

Index = ( temphum .* InnerEdge.nx  <= 0 & - temphup .* InnerEdge.nx > 0 );
fluxX( Index ) =  vfpx(Index) .* InnerEdge.nx(Index);
% wet-dry interface considered
[vfmx, vfpx, fluxX] = obj.matGetWetDryFaceVariableAndFlux( vfmx, vfpx, fluxX );

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
% [ vfpx ] = ImposeBoundaryCondition(  vfmx, vfpx, obj.EidBoundaryType);
% getCentralFluxTerm
fluxX = BoundaryEdge.nx .* ( vfmx + vfpx )./2; 
temphum = fm(:,:,2);
temphup = fp(:,:,2);

% getUpwindFluxTerm
Index = ( temphum .* BoundaryEdge.nx  > 0 & - temphup .* BoundaryEdge.nx  <= 0 );
fluxX(Index) = vfmx(Index) .* BoundaryEdge.nx(Index);

Index = ( temphum .* BoundaryEdge.nx &  - temphup .* BoundaryEdge.nx > 0 );
fluxX( Index ) =  vfpx(Index) .* BoundaryEdge.nx(Index);

% getLocalFluxTerm
fluxMx = BoundaryEdge.nx .* vfmx;

% the boundary edge contribution
UpwindedTermX = -UpwindedTermX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMx, fluxX);

[ VolumeX ] = obj.matVolumeIntegral( mesh, cell2mat(variableX));
UpwindedTermX = UpwindedTermX + VolumeX;
end

% function [ vfpx ] = ImposeBoundaryCondition(  vfmx, vfpx, EidBoundaryType)
% Index = ( EidBoundaryType == 1 );
% vfpx(Index) = vfmx(Index);
% Index = ( EidBoundaryType == -1 );
% vfpx(Index) = -vfmx(Index);
% end