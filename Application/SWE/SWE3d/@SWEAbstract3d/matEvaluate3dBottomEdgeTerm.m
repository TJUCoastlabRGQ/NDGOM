function BottomEdge_rhs = matEvaluate3dBottomEdgeTerm( obj, edge, fphys3d)
%> @brief Function used to calculate the bottom edge contribution to the RHS
%> 
%> More detailed description.
%>
%> For this case, we only consider the vertical momentum contribution to the
%> right hand side, and the contribution due to the auxialary variable is
%> not considered here. 
%>
%> @param edge The bottom edge object
%> @param fphys3d The three dimensional physical field
%>
%> @retval BottomEdge_rhs The right hand side due to the bottom edge
%> contribution

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );

%> Following two formula are used to calclulate $\mathbf F^-\cdot \mathbf n$
FluxM(:, :, 1) = fm( :, :, 1 ) .* fm( :, :, 3 ) ./ fm( :, :, 6 )  .* edge.nz;   
FluxM(:, :, 2) = fm( :, :, 2 ) .* fm( :, :, 3 ) ./ fm( :, :, 6 )  .* edge.nz;

%> Following two formula are used to calclulate $\mathbf F^+\cdot \mathbf n$
FluxP(:, :, 1) = fp( :, :, 1 ) .* fp( :, :, 3 ) ./ fp( :, :, 6 )  .* edge.nz;
FluxP(:, :, 2) = fp( :, :, 2 ) .* fp( :, :, 3 ) ./ fp( :, :, 6 )  .* edge.nz;


%> $\lambda_{max}  = max(\lambda_1,\lambda_2,\lambda_3)$, for the situation
%> here $\lambda_1 = \lambda_2 = \lambda_3$, because $n_x = n_y = 0$
lambda = abs( max( max( fm(:, :, 3) ./ fm(:, :, 6) .* edge.nz ,  - fp(:, :, 3) ./ fp(:, :, 6) .* edge.nz ) ) );

FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    bsxfun( @times,   lambda  , ( fp( :, :, 1 ) - fm( :, :, 1 ) ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    bsxfun( @times,  lambda  , ( fp( :, :, 2 ) - fm( :, :, 2 ) ) ) );

BottomEdge_rhs = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxP, FluxS );

end