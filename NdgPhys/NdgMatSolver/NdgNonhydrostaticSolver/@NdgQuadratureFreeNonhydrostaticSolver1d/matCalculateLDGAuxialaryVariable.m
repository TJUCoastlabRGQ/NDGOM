function [ qx ] = matCalculateLDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable)
 [fmx, fpx] = InnerEdge.matEvaluateSurfValue( Variable );
% Boundary condition at the wet-dry interface need to be considered
[Um, Up] = InnerEdge.matEvaluateSurfValue( Variable );
%< Inner edge contribution
fluxMX = InnerEdge.nx.*fmx; 
fluxPX = InnerEdge.nx.*fpx;
fluxSx = ( ( fmx + fpx ) ./ 2 - 1/2 .* ( Um - Up ) ) .* InnerEdge.nx;

fluxSx = obj.matGetPrimitiveVariableInnerEdgeFlux(  obj.WetDryFaceOrder, fluxSx, mesh.cell.Nfp);

termX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);

[fmx, ~] = BoundaryEdge.matEvaluateSurfValue( Variable ); 


%< Boundary edge contribution
fluxMX = BoundaryEdge.nx.*fmx;
% Boundary condition at the wet-dry interface need to be considered,
% follow matCalculateCharacteristicMatrix
fluxSX = obj.matGetPrimitiveVariableBoundaryEdgeFlux( BoundaryEdge.nx, fmx );

termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

[ VolumeIntegralX ] = obj.matVolumeIntegral( mesh, cell2mat(Variable) );

qx = termX + VolumeIntegralX;
end