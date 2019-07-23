function [qx, qy] = matCalculateLDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable)
[fmy, fpy] = InnerEdge.matEvaluateSurfValue( Variable );  [fmx, fpx] = InnerEdge.matEvaluateSurfValue( Variable );
% Boundary condition at the wet-dry interface has been considered by setting the flux at the wet-dry interface flaged by WetDryFaceOrder directly
% [Um, Up] = InnerEdge.matEvaluateSurfValue( Variable );
%< Inner edge contribution
fluxMY = InnerEdge.ny.*fmy; fluxMX = InnerEdge.nx.*fmx; 
fluxPY = InnerEdge.ny.*fpy; fluxPX = InnerEdge.nx.*fpx;
fluxSx = ( ( fmx + fpx ) ./ 2 ) .* InnerEdge.nx;
fluxSy = ( ( fmy + fpy ) ./ 2 ) .* InnerEdge.ny;
fluxSx = obj.matGetPrimitiveVariableInnerEdgeFlux(  obj.WetDryFaceOrder, fluxSx, mesh.cell.Nfp(1));
fluxSy = obj.matGetPrimitiveVariableInnerEdgeFlux(  obj.WetDryFaceOrder, fluxSy, mesh.cell.Nfp(1));

% termY = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMY, fluxPY);
% termX = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMX, fluxPX);
termY = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxPY, fluxSy);
termX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);

[fmy, ~] = BoundaryEdge.matEvaluateSurfValue( Variable );    [fmx, ~] = BoundaryEdge.matEvaluateSurfValue( Variable ); 


%< Boundary edge contribution
fluxMY = BoundaryEdge.ny.*fmy;  fluxMX = BoundaryEdge.nx.*fmx;

fluxSX = obj.matGetPrimitiveVariableBoundaryEdgeFlux( BoundaryEdge.nx, fmx );
fluxSY = obj.matGetPrimitiveVariableBoundaryEdgeFlux( BoundaryEdge.ny, fmy );

termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);
termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

[VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral( mesh, cell2mat(Variable), cell2mat(Variable));

qy = termY + VolumeIntegralY;
qx = termX + VolumeIntegralX;
end