function frhs2d_BoundarySurfaceTerm = matEvaluate2dHorizonPCEBoundaryTerm( obj, BoundaryEdge, fphys2d, fext)
%> @brief Function used to calculate the two dimentional PCE boundary surface integration term
%> 
%> More detailed description.
%>
%> @param BoundaryEdge The boundary edge object
%> @param fphys2d The two dimensional physical field
%> @param fext The exterior physical field
%>
%> @fphys2d frhs2d_BoundarySurfaceTerm Contribution to the RHS of the two dimensional PCE from the boundary face integration term


[ fm, fp ] = BoundaryEdge.matEvaluateSurfValue( fphys2d );

% apply clamped boundary condition
ind = ( BoundaryEdge.ftype == enumBoundaryCondition.Clamped );
% fp(:, ind, 4) = fext2d(:, ind, 4); 
fp(:, ind, 1) = fext(:, ind, 1);

% apply slip wall boundary condition
ind = ( BoundaryEdge.ftype == enumBoundaryCondition.SlipWall );
Hun =  fm( :, ind, 2 ) .* BoundaryEdge.nx(:, ind) + fm( :, ind, 3).* BoundaryEdge.ny(:, ind);
Hvn = -fm( :, ind, 2 ) .* BoundaryEdge.ny(:, ind) + fm( :, ind, 3).* BoundaryEdge.nx(:, ind);

fp(:, ind, 2) = - Hun .* BoundaryEdge.nx(:, ind) - Hvn .* BoundaryEdge.ny(:, ind);
fp(:, ind, 3) = - Hun .* BoundaryEdge.ny(:, ind) + Hvn .* BoundaryEdge.nx(:, ind);

lambda = max( sqrt( obj.gra .* fm(:, :, 4) ), sqrt( obj.gra .* fp(:, :, 4) ) );

FluxM = fm(:, :, 2) .* BoundaryEdge.nx + fm(:, :, 3) .* BoundaryEdge.ny;
FluxP = fp(:, :, 2) .* BoundaryEdge.nx + fp(:, :, 3) .* BoundaryEdge.ny;
FluxS = 0.5 * ( FluxM + FluxP - lambda .* ( fp(:, :, 1) - fm(:, :, 1) ) );

frhs2d_BoundarySurfaceTerm = BoundaryEdge.matEvaluateStrongFromEdgeRHS( FluxM, FluxS );
end