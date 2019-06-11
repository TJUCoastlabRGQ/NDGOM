function [ q2x ] = matCalculateLDGSecondOrderVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable, VariableX )

[fmx, fpx] = InnerEdge.matEvaluateSurfValue( VariableX ); 
[Um, Up] = InnerEdge.matEvaluateSurfValue( Variable ); 
JumpUx = InnerEdge.nx .* Um - InnerEdge.nx .* Up;
Jumpq = InnerEdge.nx .* fmx - InnerEdge.nx .* fpx;
fluxMX = InnerEdge.nx.*fmx; 
fluxPX = InnerEdge.nx.*fpx;
fluxSx = ( ( fmx + fpx ) ./ 2 - obj.IETau .* JumpUx + 1/2 .* InnerEdge.nx .* Jumpq ) .* InnerEdge.nx;

%Wet Dry Interface considered
fluxSx = obj.matGetAuxialaryVariableInnerEdgeFlux( fluxSx, Um, Up, fmx, fpx, InnerEdge.nx);

q2x = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);

[fmx, ~] = BoundaryEdge.matEvaluateSurfValue( VariableX ); 
fluxMX = BoundaryEdge.nx.*fmx;
fluxSX = obj.matGetAuxialaryVariableBoundaryEdgeFlux( Um, fmx, BoundaryEdge.nx );

q2x = - q2x - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

[VolumeIntegralX ] = obj.matVolumeIntegral( mesh, cell2mat(VariableX) );

q2x = q2x + VolumeIntegralX;
end