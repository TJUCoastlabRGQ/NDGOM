function [ W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d )
%MATEVALUATEVERTICALVELOCITY Summary of this function goes here
%   Detailed explanation goes here

% fphys3d(:, :, 6) = mesh3d.Extend2dField( fphys2d(:, :, 1) - fphys2d(:, :, 5) );

% HUt = mesh3d.Extend2dField( fphys2d(:, :, 2) ) - fphys3d(:, :, 5) .* fphys3d(:, :, 1);
% HVt = mesh3d.Extend2dField( fphys2d(:, :, 3) ) - fphys3d(:, :, 5) .* fphys3d(:, :, 2);

% dHUdx = mesh3d.rx .* ( mesh3d.cell.Dr * fphys3d(:, :, 1) ) + mesh3d.sx .* ( mesh3d.cell.Ds * fphys3d(:, :, 1) ) + mesh3d.tx .* ( mesh3d.cell.Dt * fphys3d(:, :, 1) );
% dHVdy = mesh3d.ry .* ( mesh3d.cell.Dr * fphys3d(:, :, 2) ) + mesh3d.sy .* ( mesh3d.cell.Ds * fphys3d(:, :, 2) ) + mesh3d.ty .* ( mesh3d.cell.Dt * fphys3d(:, :, 2) );

% [ dHUHV ] = obj.matEvaluateHorizontalPartialDerivativeTerm( mesh3d,data );
[ dHUHV, dHU2DHV2D ] = matEvaluateHorizontalPartialDerivativeTermNew( obj, mesh3d, fphys2d, fphys3d );
% dHVdy = obj.matEvaluatePartialDerivativeTermX( mesh, fphys3d);
% fphys2d = obj.matEvaluate2dHorizonMomentum( mesh3d, ...
%         fphys2d, fphys3d );
% fphys2d(:, :, 2) = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 1) );
% fphys2d(:, :, 3) = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 2) );

% HU2D = mesh3d.Extend2dField( fphys2d(:, :, 2) );
% HV2D = mesh3d.Extend2dField( fphys2d(:, :, 3) );

% data{1}(:,:,1) = mesh3d.Extend2dField( fphys2d(:, :, 2) );
% data{1}(:,:,2) = mesh3d.Extend2dField( fphys2d(:, :, 3) );
% [ dHU2DHV2D ] = obj.matEvaluateHorizontalPartialDerivativeTerm(mesh3d, data); 

% dHU2Ddx = mesh3d.rx .* ( mesh3d.cell.Dr * HU2D ) + mesh3d.sx .* ( mesh3d.cell.Ds * HU2D ) + mesh3d.tx .* ( mesh3d.cell.Dt * HU2D );
% dHV2Ddy = mesh3d.ry .* ( mesh3d.cell.Dr * HV2D ) + mesh3d.sy .* ( mesh3d.cell.Ds * HV2D ) + mesh3d.ty .* ( mesh3d.cell.Dt * HV2D );

W =  mesh3d.VerticalIntegralField( - dHUHV ) + dHU2DHV2D .* (1+mesh3d.z);
end

function [ Term, Term2d ] = matEvaluateHorizontalPartialDerivativeTermNew( obj, mesh3d, fphys2d, fphys3d)

    edge = mesh3d.InnerEdge;
    edge2d = obj.mesh2d(1).InnerEdge;     

    [ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
    [ fm2d, fp2d ] = edge2d.matEvaluateSurfValue( fphys2d );
    
    FluxM = fm(:,:,1) .* edge.nx + fm(:,:,2) .* edge.ny;
    FluxM2d = fm2d(:,:,2) .* edge2d.nx + fm2d(:,:,3) .* edge2d.ny;
    
    FluxP = fp(:,:,1) .* edge.nx + fp(:,:,2) .* edge.ny;
    FluxP2d = fp2d(:,:,2) .* edge2d.nx + fp2d(:,:,3) .* edge2d.ny;
        
    FluxS = obj.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4 1 2]), fp(:,:,[4 1 2]), edge );
    FluxS2d = edge.VerticalColumnIntegralField(FluxS(:,:,1));
    InnerSurface_frhs3d = edge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS(:,:,1) );
    InnerSurface_frhs2d = edge2d.matEvaluateStrongFromEdgeRHS( FluxM2d, FluxP2d, FluxS2d );
    
    edge = mesh3d.BoundaryEdge;
    edge2d = obj.mesh2d(1).BoundaryEdge;
    [ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
    [ fm2d, ~ ] = edge2d.matEvaluateSurfValue( fphys2d );

    [ fm, fp ] = obj.matImposeBoundaryCondition( edge, edge.nx, edge.ny, fm, fp, obj.fext3d{1} );

    FluxM = fm(:,:,1) .* edge.nx + fm(:,:,2) .* edge.ny;
    FluxM2d = fm2d(:,:,2) .* edge2d.nx + fm2d(:,:,3) .* edge2d.ny;
    FluxS = obj.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, fm(:,:,[4 1 2]), fp(:,:,[4 1 2]), edge );
    FluxS2d = edge.VerticalColumnIntegralField(FluxS(:,:,1));
    
    BoundarySurface_frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS(:,:,1) );
    BoundarySurface_frhs2d = edge2d.matEvaluateStrongFromEdgeRHS( FluxM2d, FluxS2d );
   
Term =...
    ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,1)) + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,1)) ) - ...
    InnerSurface_frhs3d - BoundarySurface_frhs3d  + ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,2)) + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,2)) ) ;
Term2d = obj.mesh2d.rx .* (obj.mesh2d.cell.Dr * fphys2d{1}(:,:,2)) + obj.mesh2d.sx .* (obj.mesh2d.cell.Ds * fphys2d{1}(:,:,2) ) - ...
    InnerSurface_frhs2d - BoundarySurface_frhs2d + obj.mesh2d.ry .* (obj.mesh2d.cell.Dr * fphys2d{1}(:,:,3) ) + obj.mesh2d.sy .* (obj.mesh2d.cell.Ds * fphys2d{1}(:,:,3));
Term2d = obj.meshUnion(1).Extend2dField( Term2d );
end

