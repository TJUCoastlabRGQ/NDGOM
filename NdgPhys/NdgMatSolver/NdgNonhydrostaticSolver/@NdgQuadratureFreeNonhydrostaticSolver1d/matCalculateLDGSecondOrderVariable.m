function [ q2x ] = matCalculateLDGSecondOrderVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable, VariableX )

Ux = mesh.rx .* ( mesh.cell.Dr * cell2mat(Variable) ) + mesh.sx .* ( mesh.cell.Ds * cell2mat(Variable) );
[Um, Up] = InnerEdge.matEvaluateSurfValue( Variable ); 
[Umx, Upx] = InnerEdge.matEvaluateSurfValue( num2cell(Ux, [1,2]) ); 
[fmx, fpx] = InnerEdge.matEvaluateSurfValue( VariableX ); 

JumpUx = InnerEdge.nx .* Um - InnerEdge.nx .* Up;
% Jumpq = InnerEdge.nx .* fmx - InnerEdge.nx .* fpx;
fluxMX = InnerEdge.nx.*fmx; 
fluxPX = InnerEdge.nx.*fpx;
fluxSx = ( ( Umx + Upx ) ./ 2 - obj.Tau * JumpUx ) .* InnerEdge.nx;

%Wet Dry Interface considered
fluxSx = obj.matGetAuxialaryVariableInnerEdgeFlux( fluxSx, Um, Up, Umx, Upx, InnerEdge.nx);

q2x = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);

[fmx, ~] = BoundaryEdge.matEvaluateSurfValue( VariableX ); 
[Umx, ~] = BoundaryEdge.matEvaluateSurfValue( num2cell(Ux, [1,2]) ); 
fluxMX = BoundaryEdge.nx.*fmx;
fluxSX = obj.matGetAuxialaryVariableBoundaryEdgeFlux( Um, Umx, BoundaryEdge.nx );

q2x = - q2x - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

[VolumeIntegralX ] = obj.matVolumeIntegral( mesh, cell2mat(VariableX) );

q2x = q2x + VolumeIntegralX;
end