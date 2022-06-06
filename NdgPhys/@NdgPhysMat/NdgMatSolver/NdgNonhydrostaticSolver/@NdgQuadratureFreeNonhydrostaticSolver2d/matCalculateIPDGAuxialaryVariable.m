function [qx, qy] = matCalculateIPDGAuxialaryVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable)
[fmy, fpy] = InnerEdge.matEvaluateSurfValue( Variable );  [fmx, fpx] = InnerEdge.matEvaluateSurfValue( Variable );
%< Inner edge contribution
fluxMY = InnerEdge.ny.*fmy; fluxMX = InnerEdge.nx.*fmx; 
fluxPY = InnerEdge.ny.*fpy; fluxPX = InnerEdge.nx.*fpx;
fluxSx = ( ( fmx + fpx ) ./ 2 ) .* InnerEdge.nx;
fluxSy = ( ( fmy + fpy ) ./ 2 ) .* InnerEdge.ny;

% 
faceflag = obj.matGetWetDryFace( mesh.status, InnerEdge.FToE);

WetDryFaceOrder = find( faceflag == 1);

fluxSx = obj.matGetPrimitiveVariableInnerEdgeFlux(  WetDryFaceOrder', fluxSx, mesh.cell.Nfp(1));  
fluxSy = obj.matGetPrimitiveVariableInnerEdgeFlux(  WetDryFaceOrder', fluxSy, mesh.cell.Nfp(1));

termY = InnerEdge.matEvaluateStrongFormEdgeRHS(fluxMY, fluxPY, fluxSy);
termX = InnerEdge.matEvaluateStrongFormEdgeRHS(fluxMX, fluxPX, fluxSx);

[fmy, ~] = BoundaryEdge.matEvaluateSurfValue( Variable );    
[fmx, ~] = BoundaryEdge.matEvaluateSurfValue( Variable ); 


%< Boundary edge contribution
fluxMY = BoundaryEdge.ny.*fmy;  fluxMX = BoundaryEdge.nx.*fmx;

fluxSX = obj.matGetPrimitiveVariableBoundaryEdgeFlux( BoundaryEdge.nx, fmx );  
fluxSY = obj.matGetPrimitiveVariableBoundaryEdgeFlux( BoundaryEdge.ny, fmy );

termY = - termY - BoundaryEdge.matEvaluateStrongFormEdgeRHS(fluxMY, fluxSY);
termX = - termX - BoundaryEdge.matEvaluateStrongFormEdgeRHS(fluxMX, fluxSX);

[VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral( mesh, cell2mat(Variable), cell2mat(Variable));  

qy = termY + VolumeIntegralY;
qx = termX + VolumeIntegralX;
end