function  [qx, q2x, Nonhydrop] = matAssembleCharacteristicMatrix(obj, mesh, index)
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

% [ Tempqx ]  = obj.matCalculateCharacteristicMatrix( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
% tempqy = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
[ tempqx ]  = obj.matCalculateLDGAuxialaryVariable( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]));

qx = tempqx(:);

% [ tempq2x, tempq2y ] = obj.matCalculateCharacteristicMatrix( mesh, BoundaryEdge, InnerEdge, num2cell(tempqx,[1 2]), num2cell(tempqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
% tempq2y = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(tempqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
[ tempq2x ] = obj.matCalculateLDGSecondOrderVariable( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), num2cell(tempqx,[1 2]) );

q2x = tempq2x(:);

end