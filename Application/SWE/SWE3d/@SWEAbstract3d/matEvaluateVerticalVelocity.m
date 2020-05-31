function [ W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d )
%MATEVALUATEVERTICALVELOCITY Summary of this function goes here
%   Detailed explanation goes here


    edge = mesh3d.InnerEdge;
    edge2d = obj.mesh2d(1).InnerEdge;     


    InnerSurface_frhs3d = edge.matEvaluateStrongFromEdgeRHS( obj.InnerEdgeFluxM3d{1}(:,:,1), obj.InnerEdgeFluxP3d{1}(:,:,1), obj.InnerEdgeFluxS3d{1}(:,:,1) );
    InnerSurface_frhs2d = edge2d.matEvaluateStrongFromEdgeRHS( obj.InnerEdgeFluxM2d{1}(:,:,1), obj.InnerEdgeFluxP2d{1}(:,:,1), obj.InnerEdgeFluxS2d{1}(:,:,1) );
    
    edge = mesh3d.BoundaryEdge;
    edge2d = obj.mesh2d(1).BoundaryEdge;
    
    BoundarySurface_frhs3d = edge.matEvaluateStrongFormEdgeRHS( obj.BoundaryEdgeFluxM3d{1}(:,:,1), obj.BoundaryEdgeFluxS3d{1}(:,:,1) );
    BoundarySurface_frhs2d = edge2d.matEvaluateStrongFromEdgeRHS( obj.BoundaryEdgeFluxM2d{1}(:,:,1), obj.BoundaryEdgeFluxS2d{1}(:,:,1) );
   
dHUHV =...
    ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,1)) + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,1)) ) - ...
    InnerSurface_frhs3d - BoundarySurface_frhs3d  + ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,2)) + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,2)) ) ;
% Term2d = ;
dHU2DHV2D = obj.meshUnion(1).Extend2dField( obj.mesh2d.rx .* (obj.mesh2d.cell.Dr * fphys2d{1}(:,:,2)) + obj.mesh2d.sx .* (obj.mesh2d.cell.Ds * fphys2d{1}(:,:,2) ) - ...
    InnerSurface_frhs2d - BoundarySurface_frhs2d + obj.mesh2d.ry .* (obj.mesh2d.cell.Dr * fphys2d{1}(:,:,3) ) + obj.mesh2d.sy .* (obj.mesh2d.cell.Ds * fphys2d{1}(:,:,3)) );

W =  mesh3d.VerticalIntegralField( -1 * dHUHV ) + dHU2DHV2D .* (1+mesh3d.z);
end

