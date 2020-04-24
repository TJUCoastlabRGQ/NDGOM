function [ Term ] = matEvaluateHorizontalPartialDerivativeTerm(obj, mesh3d, fphys3d)

    edge = mesh3d.InnerEdge;

    [ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
    
    FluxM = fm(:,:,1) .* edge.nx + fm(:,:,2) .* edge.ny;
    
    FluxP = fp(:,:,1) .* edge.nx + fp(:,:,2) .* edge.ny;
    
    FluxS = obj.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4 1 2]), fp(:,:,[4 1 2]), edge );
    
    InnerSurface_frhs3d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS(:,:,1) );
        
    edge = mesh3d.BoundaryEdge;
    [ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );

    ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
    Hun =  fm( :, ind, 1 ) .* edge.nx(:, ind) + fm( :, ind, 2).* edge.ny(:, ind);
    Hvn = -fm( :, ind, 1 ) .* edge.ny(:, ind) + fm( :, ind, 2).* edge.nx(:, ind);

    fp(:, ind, 1) = - Hun .* edge.nx(:, ind) - Hvn .* edge.ny(:, ind);
    fp(:, ind, 2) = - Hun .* edge.ny(:, ind) + Hvn .* edge.nx(:, ind);

    FluxM = fm(:,:,1) .* edge.nx + fm(:,:,2) .* edge.ny;

    FluxS = obj.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4 1 2]), fp(:,:,[4 1 2]), edge );

    
    BoundarySurface_frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS(:,:,1) );

Term =...
    ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,1)) + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,1)) + mesh3d.tx .* (mesh3d.cell.Dt * fphys3d{1}(:,:,1)) - ...
    InnerSurface_frhs3d - BoundarySurface_frhs3d ) + ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,2)) + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,2)) + ...
    mesh3d.ty .* (mesh3d.cell.Dt * fphys3d{1}(:,:,2)) ) ;

end