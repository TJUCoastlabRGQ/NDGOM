function [ qx ] = matCalculateIPDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable)
[fmx, fpx] = InnerEdge.matEvaluateSurfValue( Variable );
%< Inner edge contribution
 fluxMX = InnerEdge.nx.*fmx; 
 fluxPX = InnerEdge.nx.*fpx;
fluxSx = ( ( fmx + fpx ) ./ 2 ) .* InnerEdge.nx;

% 
faceflag = obj.matGetWetDryFace( mesh.status, InnerEdge.FToE);

WetDryFaceOrder = find( faceflag == 1);

fluxSx = obj.matGetPrimitiveVariableInnerEdgeFlux(  WetDryFaceOrder', fluxSx, mesh.cell.Nfp(1));  

termX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);

[fmx, ~] = BoundaryEdge.matEvaluateSurfValue( Variable ); 


%< Boundary edge contribution
 fluxMX = BoundaryEdge.nx.*fmx;

fluxSX = obj.matGetPrimitiveVariableBoundaryEdgeFlux( BoundaryEdge.nx, fmx );  

termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

[VolumeIntegralX] = obj.matVolumeIntegral( mesh, cell2mat(Variable) );  

qx = termX + VolumeIntegralX;
end