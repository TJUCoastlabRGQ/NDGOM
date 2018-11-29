function matEvaluateRHS( obj, fphys2d, fphys3d )
%MATEVALUATERHS Summary of this function goes here
%   Detailed explanation goes here

for m = 1:obj.Nmesh
    mesh2d = obj.mesh2d(m);
    mesh3d = obj.mesh3d(m);
    
    % evaluate 2d PCE volume integral term
    fphys2d{m} = obj.matEvaluate2dHorizonMomentum( mesh3d, ...
        fphys2d{m}, fphys3d{m} );
    
    obj.frhs2d{m} = obj.matEvaluate2dHorizonPCEVolumeTerm( mesh2d, fphys2d{m} );

    % evaluate 2d PCE surface integral term
    obj.frhs2d{m} = obj.frhs2d{m} + obj.matEvaluate2dHorizonPCEInnerSurfaceTerm( ...
        mesh2d.InnerEdge, fphys2d );
    
    obj.frhs2d{m} = obj.frhs2d{m} + obj.matEvaluate2dHorizonPCEBoundaryTerm( ...
        obj, mesh2d.BoundaryEdge, fphys2d, obj.fext2d{m});
    
    % extend 3d field
    fphys3d = obj.matEvaluate3dAuxiliaryVariable( obj, mesh3d, fphys3d );   %1
    
    % evaluate 3d velocity volume integral
    obj.frhs3d{m} = obj.matEvaluate3dVolumeTerm( mesh3d, fphys3d{m} );   %2
    
    % evaluate 3d velocity horizontal surface integral
    obj.frhs3d{m} = obj.frhs3d{m} +obj.matEvaluate3dInnerSurfaceTerm( ...           %3
         mesh3d.InnerEdge, fphys3d );
    
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_HorizontalBoundaryKernel( ...     %4
        obj, mesh3d.BoundaryEdge, fphys3d, obj.fext3d{m} );
    
    % evaluate 3d velocity field surface integral
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_SurfaceBoundaryKernel( ...        %5
        obj, mesh3d.SurfaceBoundaryEdge, fphys3d );
    
    % evaluate 3d velocity field bottom integral
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_BottomSurfaceKernel( ...         %6
        obj, mesh3d.BottomEdge, fphys3d);
    
    obj.frhs3d{m} = obj.frhs3d{m} + Evaluate3d_BottomBoundaryKernel( ...       %7
        obj, mesh3d.BottomBoundaryEdge, fphys3d );
end

end

function frhs3d = Evaluate3d_SideSurfaceKernel( gra, edge, fphys3d )

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = gra .* fm( :, :, 7 ) .* edge.nx;
FluxP(:, :, 1) = gra .* fp( :, :, 7 ) .* edge.nx;
FluxM(:, :, 2) = gra .* fm( :, :, 7 ) .* edge.ny;
FluxP(:, :, 2) = gra .* fp( :, :, 7 ) .* edge.ny;

lambda = max( sqrt( gra .* fm(:, :, 6) ), sqrt( gra .* fp(:, :, 6) ) );

FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    lambda .* ( fp( :, :, 1 ) - fm( :, :, 1 ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    lambda .* ( fp( :, :, 2 ) - fm( :, :, 2 ) ) );

frhs3d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS );
end

function frhs3d = Evaluate3d_SurfaceBoundaryKernel( obj, edge, fphys3d )

Hmiu = sqrt( obj.miu0 );

[ fm, ~ ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = Hmiu .* fm( :, :, 4 ) .* edge.nz;
FluxM(:, :, 2) = Hmiu .* fm( :, :, 5 ) .* edge.nz;

FluxS(:, :, 1) = zeros(size( fm( :, :, 4 )));
FluxS(:, :, 2) = zeros(size( fm( :, :, 5 )));


frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );
end

function frhs3d = Evaluate3d_BottomSurfaceKernel( ...
    obj, edge, fphys3d )

Hmiu = sqrt( obj.miu0 );
tau = 1e-3;

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = Hmiu .* ( fm( :, :, 4 ) .* edge.nz + tau * fm(:, :, 1) );
FluxM(:, :, 2) = Hmiu .* ( fm( :, :, 5 ) .* edge.nz + tau * fm(:, :, 2) );
FluxP(:, :, 1) = Hmiu .* ( fp( :, :, 4 ) .* edge.nz + tau * fp(:, :, 1) );
FluxP(:, :, 2) = Hmiu .* ( fp( :, :, 5 ) .* edge.nz + tau * fp(:, :, 2) );

frhs3d = edge.matEvaluateStrongFormEdgeCentralRHS( FluxM, FluxP );
end

function frhs3d = Evaluate3d_BottomBoundaryKernel( obj, edge, fphys3d )
Hmiu = sqrt( obj.miu0 );

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = Hmiu .* fm( :, :, 4 ) .* edge.nz;
FluxM(:, :, 2) = Hmiu .* fm( :, :, 5 ) .* edge.nz;

FluxS(:, :, 1) =  - obj.K .* fp( :, :, 1 ) ./ fp( :, :, 6 ) .* edge.nz;
FluxS(:, :, 2) =  - obj.K .* fp( :, :, 2 ) ./ fp( :, :, 6 ) .* edge.nz;

frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );
end
