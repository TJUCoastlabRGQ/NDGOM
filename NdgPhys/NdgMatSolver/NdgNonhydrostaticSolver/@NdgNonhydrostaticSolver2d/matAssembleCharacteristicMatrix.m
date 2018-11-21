function  [qx, qy, q2x, q2y, qbx, qby, fqbx, fqby] = matAssembleCharacteristicMatrix(obj, mesh, index)
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
%< @param[out] qbx product of nonhydrostatic pressure with topography gradient bx in x direction
%< @param[out] qby product of nonhydrostatic pressure with topography gradient by in y direction
%< @param[out] fqbx product of nonhydrostatic pressure with topography gradient bx in x direction, this term is further included in the flux term
%< @param[out] fqby product of nonhydrostatic pressure with topography gradient by in y direction, this term is further included in the flux term
K = mesh.K; Np = mesh.cell.Np;
Nonhydro = zeros( Np, K );
gmat = zeros( Np, K );
gmat(index) = 1.0/2.0;
Nonhydro(index) = 1;
qbx = zeros(K*Np,1);qby = zeros(K*Np,1);

InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
tempqx = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
tempqy = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(gmat,[1 2]), enumNonhydroBoundaryCondition.Zero);
qx = tempqx(:);
qy = tempqy(:);
tempq2x = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(tempqx,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
tempq2y = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(tempqy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);
q2x = tempq2x(:); q2y = tempq2y(:);

tempfqbx = obj.matCalculateCharacteristicMatrixX( mesh, BoundaryEdge, InnerEdge, num2cell(Nonhydro.*obj.bx,[1 2]), enumNonhydroBoundaryCondition.Zero);
tempfqby = obj.matCalculateCharacteristicMatrixY( mesh, BoundaryEdge, InnerEdge, num2cell(Nonhydro.*obj.by,[1 2]), enumNonhydroBoundaryCondition.Zero);
fqbx = tempfqbx(:);
fqby = tempfqby(:);
qbx(index) = obj.bx(index);qby(index) = obj.by(index);
end