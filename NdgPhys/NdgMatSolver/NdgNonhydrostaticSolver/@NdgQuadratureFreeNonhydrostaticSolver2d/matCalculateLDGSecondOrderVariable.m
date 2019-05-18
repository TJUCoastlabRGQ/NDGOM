function [q2x, q2y] = matCalculateLDGSecondOrderVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable, VariableX, VariableY )
% C11 = 100;
[fmy, fpy] = InnerEdge.matEvaluateSurfValue( VariableY );  
[fmx, fpx] = InnerEdge.matEvaluateSurfValue( VariableX ); 
[Um, Up] = InnerEdge.matEvaluateSurfValue( Variable ); 
JumpUx = InnerEdge.nx .* Um - InnerEdge.nx .* Up;
JumpUy = InnerEdge.ny .* Um - InnerEdge.ny .* Up;
Jumpq = InnerEdge.nx .* fmx + InnerEdge.ny .* fmy...
    - InnerEdge.nx .* fpx - InnerEdge.ny .* fpy;
fluxMY = InnerEdge.ny.*fmy; fluxMX = InnerEdge.nx.*fmx; 
fluxPY = InnerEdge.ny.*fpy; fluxPX = InnerEdge.nx.*fpx;
fluxSx = ( ( fmx + fpx ) ./ 2 - obj.IETau .* JumpUx + 1/2 .* InnerEdge.nx .* Jumpq ) .* InnerEdge.nx;
fluxSy = ( ( fmy + fpy ) ./ 2 - obj.IETau .* JumpUy + 1/2 .* InnerEdge.ny .* Jumpq ) .* InnerEdge.ny;

q2x = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);
q2y = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxPY, fluxSy);

[fmy, ~] = BoundaryEdge.matEvaluateSurfValue( VariableY );    [fmx, ~] = BoundaryEdge.matEvaluateSurfValue( VariableX ); 
fluxMY = BoundaryEdge.ny.*fmy;  fluxMX = BoundaryEdge.nx.*fmx;
fluxSX = zeros( size(fluxMX) );
fluxSY = zeros( size(fluxMY) );

q2x = - q2x - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
q2y = - q2y - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);

[VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral( mesh, cell2mat(VariableX), cell2mat(VariableY) );

q2y = q2y + VolumeIntegralY;
q2x = q2x + VolumeIntegralX;
end