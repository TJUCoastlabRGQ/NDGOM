function BottomBoundary_rhs = matEvaluate3dBottomBoundaryTerm( obj, edge, fphys3d)
%> @brief Function used to calculate the bottom boundary contribution to the RHS
%> 
%> More detailed description.
%>
%> For this case, we only consider the vertical momentum contribution to the
%> right hand side, and the contribution due to the auxialary variable is
%> not considered here. For this considered part, the zero boundary
%> condition is considered. i.e. '$\omega^+ = -\omega^-$'. We do nothing
%> here as the Nuemann boundary condition is applied for the primitive
%> variable '$hu$' and '$hv$', while the homogeneous Dirichlet boundary condition
%> for '$\omega$'. This implies the Local Lax-Friencher flux is zero.
%>
%> @param edge The bottom boundary edge object
%> @param fphys3d The three dimensional physical field
%>
%> @retval BottomBoundary_rhs The right hand side due to the bottom boundary contribution


[ fm, ~ ] = edge.matEvaluateSurfValue( fphys3d );

FluxM(:, :, 1) = fm(:,:,1) .* fm(:,:,3) ./ fm(:,:,6) .* edge.nz;
FluxM(:, :, 2) = fm(:,:,2) .* fm(:,:,3) ./ fm(:,:,6) .* edge.nz;

FluxS(:, :, 1) = zeros(size( fm( :, :, 1 )));
FluxS(:, :, 2) = zeros(size( fm( :, :, 1 )));

BottomBoundary_rhs = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );

end