function [ Term ] = matEvaluateHorizontalPartialDerivativeTerm(obj, mesh3d, fphys3d)

    edge = mesh3d.InnerEdge;
    edge2d = obj.mesh2d(1).InnerEdge;

    [ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
    
    FluxM = fm(:,:,1) .* edge.nx + fm(:,:,2) .* edge.ny;
    FluxM2d = edge.VerticalColumnIntegralField(FluxM);
    
    FluxP = fp(:,:,1) .* edge.nx + fp(:,:,2) .* edge.ny;
    FluxP2d = edge.VerticalColumnIntegralField(FluxP);
    
    FluxS = obj.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4 1 2]), fp(:,:,[4 1 2]), edge );
    FluxS2d = edge.VerticalColumnIntegralField(FluxS(:,:,1));
    InnerSurface_frhs3d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS(:,:,1) );
    InnerSurface_frhs2d = edge2d.matEvaluateStrongFromEdgeRHS( FluxM2d, FluxP2d, FluxS2d );
    
    edge = mesh3d.BoundaryEdge;
    [ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );

    ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
    Hun =  fm( :, ind, 1 ) .* edge.nx(:, ind) + fm( :, ind, 2).* edge.ny(:, ind);
    Hvn = -fm( :, ind, 1 ) .* edge.ny(:, ind) + fm( :, ind, 2).* edge.nx(:, ind);

    fp(:, ind, 1) = - Hun .* edge.nx(:, ind) - Hvn .* edge.ny(:, ind);
    fp(:, ind, 2) = - Hun .* edge.ny(:, ind) + Hvn .* edge.nx(:, ind);

    FluxM = fm(:,:,1) .* edge.nx + fm(:,:,2) .* edge.ny;
    FluxM2d = edge.VerticalColumnIntegralField(FluxM);
    FluxS = obj.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4 1 2]), fp(:,:,[4 1 2]), edge );
    FluxS2d = edge.VerticalColumnIntegralField(FluxS(:,:,1));
    
    BoundarySurface_frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS(:,:,1) );
    BoundarySurface_frhs2d = edge.matEvaluateStrongFormEdgeRHS( FluxM2d, FluxS2d );
    
    [ HU ] = mesh.VerticalColumnIntegralField( fphys3d{m}(:, :, 1) );
    [ HV ] = mesh.VerticalColumnIntegralField( fphys3d{m}(:, :, 2) );    
Term =...
    ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,1)) + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,1)) ) - ...
    InnerSurface_frhs3d - BoundarySurface_frhs3d  + ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,2)) + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,2)) ) ;
Term2d = mesh2d.rx .* (mesh2d.cell.Dr * HU) + mesh2d.sx .* (mesh2d.cell.Ds * HU) - ...
    InnerSurface_frhs2d - BoundarySurface_frhs2d + mesh2d.ry .* (mesh2d.cell.Dr * HV) + mesh2d.sy .* (mesh2d.cell.Ds * HV);
end