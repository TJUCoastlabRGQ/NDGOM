function [q2x, q2y] = matCalculateLDGSecondOrderVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable, VariableX, VariableY )
% C11 = 100;
Ux = mesh.rx .* ( mesh.cell.Dr * cell2mat(Variable) ) + mesh.sx .* ( mesh.cell.Ds * cell2mat(Variable) );
Uy = mesh.ry .* ( mesh.cell.Dr * cell2mat(Variable) ) + mesh.sy .* ( mesh.cell.Ds * cell2mat(Variable) );
[ Umx, Upx ] = InnerEdge.matEvaluateSurfValue( num2cell(Ux,[1,2]) );  
[ Umy, Upy ] = InnerEdge.matEvaluateSurfValue( num2cell(Uy,[1,2]) );  

[fmy, fpy] = InnerEdge.matEvaluateSurfValue( VariableY );  
[fmx, fpx] = InnerEdge.matEvaluateSurfValue( VariableX ); 
[Um, Up] = InnerEdge.matEvaluateSurfValue( Variable ); 
JumpUx = InnerEdge.nx .* Um - InnerEdge.nx .* Up;
JumpUy = InnerEdge.ny .* Um - InnerEdge.ny .* Up;
% Jumpq = InnerEdge.nx .* fmx + InnerEdge.ny .* fmy...
%     - InnerEdge.nx .* fpx - InnerEdge.ny .* fpy;
fluxMY = InnerEdge.ny.*fmy; fluxMX = InnerEdge.nx.*fmx; 
fluxPY = InnerEdge.ny.*fpy; fluxPX = InnerEdge.nx.*fpx;

fluxSx = ( ( Umx + Upx )./2 - obj.Tau * JumpUx ) .* InnerEdge.nx;
fluxSy = ( ( Umy + Upy )./2 - obj.Tau * JumpUy ) .* InnerEdge.ny;

%Wet Dry Interface considered
fluxSx = obj.matGetAuxialaryVariableInnerEdgeFlux( fluxSx, Um, Up, Umx, Upx, InnerEdge.nx);
fluxSy = obj.matGetAuxialaryVariableInnerEdgeFlux( fluxSy, Um, Up, Umy, Upy, InnerEdge.ny);

q2x = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);
q2y = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxPY, fluxSy);

[fmy, ~] = BoundaryEdge.matEvaluateSurfValue( VariableY );    [fmx, ~] = BoundaryEdge.matEvaluateSurfValue( VariableX ); 
fluxMY = BoundaryEdge.ny.*fmy;  fluxMX = BoundaryEdge.nx.*fmx;
[Umx, ~] = BoundaryEdge.matEvaluateSurfValue( num2cell(Ux, [1,2]) ); 
[Umy, ~] = BoundaryEdge.matEvaluateSurfValue( num2cell(Uy, [1,2]) ); 
fluxSX = obj.matGetAuxialaryVariableBoundaryEdgeFlux( Um, Umx, BoundaryEdge.nx );
fluxSY = obj.matGetAuxialaryVariableBoundaryEdgeFlux( Um, Umy, BoundaryEdge.ny );

q2x = - q2x - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
q2y = - q2y - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);

[VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral( mesh, cell2mat(VariableX), cell2mat(VariableY) );

q2y = q2y + VolumeIntegralY;
q2x = q2x + VolumeIntegralX;
end