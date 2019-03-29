function [ TermX, TermY ] = matEvaluatePartialDerivativeTermX(obj, mesh3d, fphys3d)

    edge = mesh3d.InnerEdge;

    [ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
    
    FluxMx = fm(:,:,1) .* edge.nx; FluxMy = fm(:,:,2) .* edge.ny;
    
    FluxPx = fp(:,:,1) .* edge.nx; FluxPy = fp(:,:,2) .* edge.ny;
    
    FluxSx = 0.5 .* ( FluxMx + FluxPx );
    FluxSy = 0.5 .* ( FluxMy + FluxPy );

    InnerSurface_frhs3d_x = edge.matEvaluateStrongFromEdgeRHS( FluxMx, FluxPx, FluxSx );
    InnerSurface_frhs3d_y = edge.matEvaluateStrongFromEdgeRHS( FluxMy, FluxPy, FluxSy );
    
    clear FluxMx FluxMy FluxPx FluxPy;
    
    edge = mesh3d.BoundaryEdge;
    [ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );

    ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
    Hun =  fm( :, ind, 1 ) .* edge.nx(:, ind) + fm( :, ind, 2).* edge.ny(:, ind);
    Hvn = -fm( :, ind, 1 ) .* edge.ny(:, ind) + fm( :, ind, 2).* edge.nx(:, ind);

    fp(:, ind, 1) = - Hun .* edge.nx(:, ind) - Hvn .* edge.ny(:, ind);
    fp(:, ind, 2) = - Hun .* edge.ny(:, ind) + Hvn .* edge.nx(:, ind);

    FluxMx = fm(:,:,1) .* edge.nx; FluxMy = fm(:,:,2) .* edge.ny;  
    FluxPx = fp(:,:,1) .* edge.nx; FluxPy = fp(:,:,2) .* edge.ny;    

    FluxSx = 0.5 .* ( FluxMx + FluxPx );
    FluxSy = 0.5 .* ( FluxMy + FluxPy );    
    
    BoundarySurface_frhs3d_x = edge.matEvaluateStrongFormEdgeRHS( FluxMx, FluxSx );
    BoundarySurface_frhs3d_y = edge.matEvaluateStrongFormEdgeRHS( FluxMy, FluxSy );

TermX =...
    ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,1)) + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,1)) + mesh3d.tx .* (mesh3d.cell.Dt * fphys3d{1}(:,:,1)) - ...
    InnerSurface_frhs3d_x - BoundarySurface_frhs3d_x ) ;

TermY =...
     ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,2)) + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,2)) + mesh3d.ty .* (mesh3d.cell.Dt * fphys3d{1}(:,:,2)) - ...
    InnerSurface_frhs3d_y - BoundarySurface_frhs3d_y ) ;
end