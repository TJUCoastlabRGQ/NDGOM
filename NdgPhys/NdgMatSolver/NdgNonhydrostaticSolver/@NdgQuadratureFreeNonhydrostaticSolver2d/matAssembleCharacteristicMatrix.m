function  [qx, qy, q2x, q2y, Nonhydrop] = matAssembleCharacteristicMatrix(obj, mesh, index)
%<@brief function used to calculate the term used when assemble the global stiff matrix
%<@detail In this version, the zero boundary condition for nonhydrostatic pressure is imposed at
%< the general clamped boundary, and the zero grad boundary condition at the wall boundary.
%< Also the zero gradient boundary condition for the auxialary variable is imposed at the
%< general clamped boundary, and the zero boundary condition at the wall boundary
%< @param[in] mesh The mesh object
%< @param[in] index The number of the studied point
%< @param[out] qx first derivative of nonhydrostatic pressure with respect to x
%< @param[out] qy first derivative of nonhydrostatic pressure with respect to y
%< @param[out] q2x second derivative of nonhydrostatic pressure with respect to x
%< @param[out] q2y second derivative of nonhydrostatic pressure with respect to y

K = mesh.K; Np = mesh.cell.Np;
Nonhydro = zeros( Np, K );
gmat = zeros( Np, K );
gmat(index) = 1;
Nonhydro(index) = 1;
Nonhydrop  =  Nonhydro(:);

InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;

[ Tempqx, Tempqy ]  = obj.matCalculateCharacteristicMatrix( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
% tempqy = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
[ tempqx, tempqy ]  = obj.matCalculateLDGAuxialaryVariable( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]));

qx = Tempqx(:);
qy = Tempqy(:);

% [ tempq2x, tempq2y ] = obj.matCalculateCharacteristicMatrix( mesh, BoundaryEdge, InnerEdge, num2cell(tempqx,[1 2]), num2cell(tempqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
% tempq2y = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(tempqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
[ tempq2x, tempq2y ] = obj.matCalculateLDGSecondOrderVariable( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), num2cell(tempqx,[1 2]), num2cell(tempqy,[1 2]) );

q2x = tempq2x(:); q2y = tempq2y(:);

end

% function [q2x, q2y] = getSecondOrderTerm( obj, mesh, BoundaryEdge, InnerEdge, Variable, VariableX, VariableY )
% C11 = 1;
% [fmy, fpy] = InnerEdge.matEvaluateSurfValue( VariableY );  
% [fmx, fpx] = InnerEdge.matEvaluateSurfValue( VariableX ); 
% [Um, Up] = InnerEdge.matEvaluateSurfValue( Variable ); 
% JumpUx = InnerEdge.nx .* Um - InnerEdge.nx .* Up;
% JumpUy = InnerEdge.ny .* Um - InnerEdge.ny .* Up;
% Jumpq = InnerEdge.nx .* fmx + InnerEdge.ny .* fmy...
%     - InnerEdge.nx .* fpx - InnerEdge.ny .* fpy;
% fluxMY = InnerEdge.ny.*fmy; fluxMX = InnerEdge.nx.*fmx; 
% fluxPY = InnerEdge.ny.*fpy; fluxPX = InnerEdge.nx.*fpx;
% fluxSx = ( ( fmx + fpx ) ./ 2 - C11 .* JumpUx + 1/2 .* InnerEdge.nx .* Jumpq ) .* InnerEdge.nx;
% fluxSy = ( ( fmy + fpy ) ./ 2 - C11 .* JumpUy + 1/2 .* InnerEdge.ny .* Jumpq ) .* InnerEdge.ny;
% 
% q2x = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxPY, fluxSy);
% q2y = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);
% 
% [fmy, ~] = BoundaryEdge.matEvaluateSurfValue( VariableY );    [fmx, ~] = BoundaryEdge.matEvaluateSurfValue( VariableX ); 
% fluxMY = BoundaryEdge.ny.*fmy;  fluxMX = BoundaryEdge.nx.*fmx;
% fluxSX = zeros( size(fluxMX) );
% fluxSY = zeros( size(fluxMY) );
% 
% q2x = - q2x - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
% q2y = - q2y - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);
% 
% [VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral( mesh, cell2mat(VariableX), cell2mat(VariableY) );
% 
% q2y = q2y + VolumeIntegralY;
% q2x = q2x + VolumeIntegralX;
% end

% function [qx, qy] = getAuxialiaryVariable( obj, mesh, BoundaryEdge, InnerEdge, Variable)
% [fmy, fpy] = InnerEdge.matEvaluateSurfValue( Variable );  [fmx, fpx] = InnerEdge.matEvaluateSurfValue( Variable );   
% [Um, Up] = InnerEdge.matEvaluateSurfValue( Variable );
% %< Inner edge contribution
% fluxMY = InnerEdge.ny.*fmy; fluxMX = InnerEdge.nx.*fmx; 
% fluxPY = InnerEdge.ny.*fpy; fluxPX = InnerEdge.nx.*fpx;
% fluxSx = ( ( fmx + fpx ) ./ 2 - 1/2 .* ( Um - Up ) ) .* InnerEdge.nx;
% fluxSy = ( ( fmy + fpy ) ./ 2 - 1/2 .* ( Um - Up ) ) .* InnerEdge.ny;
% % termY = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMY, fluxPY);
% % termX = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMX, fluxPX);
% termY = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxPY, fluxSy);
% termX = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxPX, fluxSx);
% 
% [fmy, ~] = BoundaryEdge.matEvaluateSurfValue( Variable );    [fmx, ~] = BoundaryEdge.matEvaluateSurfValue( Variable ); 
% 
% 
% %< Boundary edge contribution
% fluxMY = BoundaryEdge.ny.*fmy;  fluxMX = BoundaryEdge.nx.*fmx;
% % fluxSY = BoundaryEdge.ny.*(fpy + fmy)./2; fluxSX = BoundaryEdge.nx.*(fpx + fmx)./2;
% 
% fluxSX = ( fmx ) .* BoundaryEdge.nx;
% fluxSY = ( fmy ) .* BoundaryEdge.ny;
% 
% termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);
% termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
% 
% [VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral( mesh, cell2mat(Variable), cell2mat(Variable));
% 
% qy = termY + VolumeIntegralY;
% qx = termX + VolumeIntegralX;
% end


