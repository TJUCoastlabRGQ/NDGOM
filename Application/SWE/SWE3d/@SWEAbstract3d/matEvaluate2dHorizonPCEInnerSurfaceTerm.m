function frhs2d_InnerSurfaceTerm = matEvaluate2dHorizonPCEInnerSurfaceTerm( obj, gra, InnerEdge, fphys2d)


[ fm, fp ] = InnerEdge.matEvaluateSurfValue( fphys2d );
lambda = max( sqrt( gra .* fm(:, :, 4) ), sqrt( gra .* fp(:, :, 4) ) );

FluxM = fm(:, :, 2) .* InnerEdge.nx + fm(:, :, 3) .* InnerEdge.ny;
FluxP = fp(:, :, 2) .* InnerEdge.nx + fp(:, :, 3) .* InnerEdge.ny;
FluxS = 0.5 * ( FluxM + FluxP - lambda .* ( fp(:, :, 1) - fm(:, :, 1)  ) );

frhs2d_InnerSurfaceTerm = InnerEdge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS );

end