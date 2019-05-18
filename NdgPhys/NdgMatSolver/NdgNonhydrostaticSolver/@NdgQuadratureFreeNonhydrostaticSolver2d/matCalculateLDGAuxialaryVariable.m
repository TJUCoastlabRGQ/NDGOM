function [qx, qy] = matCalculateLDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable)
[fmy, fpy] = InnerEdge.matEvaluateSurfValue( Variable );  [fmx, fpx] = InnerEdge.matEvaluateSurfValue( Variable );
% Boundary condition at the wet-dry interface need to be considered
[Um, Up] = InnerEdge.matEvaluateSurfValue( Variable );
%< Inner edge contribution
fluxMY = InnerEdge.ny.*fmy; fluxMX = InnerEdge.nx.*fmx; 
fluxPY = InnerEdge.ny.*fpy; fluxPX = InnerEdge.nx.*fpx;
fluxSx = ( ( fmx + fpx ) ./ 2 - 1/2 .* ( Um - Up ) ) .* InnerEdge.nx;
fluxSy = ( ( fmy + fpy ) ./ 2 - 1/2 .* ( Um - Up ) ) .* InnerEdge.ny;
% termY = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMY, fluxPY);
% termX = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMX, fluxPX);
termY = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxPY, fluxSy);
termX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);

[fmy, ~] = BoundaryEdge.matEvaluateSurfValue( Variable );    [fmx, ~] = BoundaryEdge.matEvaluateSurfValue( Variable ); 


%< Boundary edge contribution
fluxMY = BoundaryEdge.ny.*fmy;  fluxMX = BoundaryEdge.nx.*fmx;
% Boundary condition at the wet-dry interface need to be considered,
% follow matCalculateCharacteristicMatrix
fluxSX = ( fmx ) .* BoundaryEdge.nx;
fluxSY = ( fmy ) .* BoundaryEdge.ny;

termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);
termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

[VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral( mesh, cell2mat(Variable), cell2mat(Variable));

qy = termY + VolumeIntegralY;
qx = termX + VolumeIntegralX;
end