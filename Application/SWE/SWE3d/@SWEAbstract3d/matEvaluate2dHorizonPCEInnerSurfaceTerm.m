function frhs2d_InnerSurfaceTerm = matEvaluate2dHorizonPCEInnerSurfaceTerm( obj, InnerEdge, fphys2d)
%> @brief Function used to calculate the two dimentional PCE inner surface term
%> 
%> More detailed description.
%>
%> @param InnerEdge The inner edge object
%> @param fphys2d The two dimensional physical field
%>
%> @fphys2d frhs2d_InnerSurfaceTerm Contribution to the RHS of the two dimensional PCE from the inner face integration term

[ fm, fp ] = InnerEdge.matEvaluateSurfValue( fphys2d );
%> $\lambda = abs( max(sqrt{(gH^_)},sqrt{(gH^+)}))$
lambda = abs( max( max( sqrt( obj.gra .* fm(:, :, 4) ), sqrt( obj.gra .* fp(:, :, 4) ) )) );

FluxM = fm(:, :, 2) .* InnerEdge.nx + fm(:, :, 3) .* InnerEdge.ny;
FluxP = fp(:, :, 2) .* InnerEdge.nx + fp(:, :, 3) .* InnerEdge.ny;
%> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
FluxS = 0.5 * ( FluxM + FluxP - bsxfun( @times, lambda, ( fp(:, :, 4) - fm(:, :, 4)  ) ) );

frhs2d_InnerSurfaceTerm = InnerEdge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS );

end