function [ Omega , W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d, time )
%MATEVALUATEVERTICALVELOCITY Summary of this function goes here
%   Detailed explanation goes here
% [~,~,~,Omega,~] = obj.matGetExactSolution( mesh3d.x, mesh3d.y, mesh3d.z, time);
% W = zeros(size(Omega));
edge = mesh3d.InnerEdge;
edge2d = obj.mesh2d(1).InnerEdge;
[ InnerEdgefm3d{1}, InnerEdgefp3d{1} ] = edge.matEvaluateSurfValue( fphys3d );
[ InnerEdgefm2d{1}, InnerEdgefp2d{1} ] = edge2d.matEvaluateSurfValue( fphys2d );
[ InnerEdgeFluxM3d{1} ] = obj.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, InnerEdgefm3d{1} );
[ InnerEdgeFluxM2d{1} ] = InnerEdgefm2d{1}(:,:,2) .* edge2d.nx + InnerEdgefm2d{1}(:,:,3) .* edge2d.ny;
[ InnerEdgeFluxP3d{1} ] = obj.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, InnerEdgefp3d{1} );
[ InnerEdgeFluxP2d{1} ] = InnerEdgefp2d{1}(:,:,2) .* edge2d.nx + InnerEdgefp2d{1}(:,:,3) .* edge2d.ny;
[ InnerEdgeFluxS3d{1} ] = obj.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, InnerEdgefm3d{1}(:,:,[4, 1, 2]), InnerEdgefp3d{1}(:,:,[4, 1, 2]), edge );
[ InnerEdgeFluxS2d{1}(:,:,1) ] = edge.VerticalColumnIntegralField(InnerEdgeFluxS3d{1}(:,:,1) );

[ InnerSurface_frhs3d ] = edge.matEvaluateStrongFromEdgeRHS( InnerEdgeFluxM3d{1}(:,:,1), InnerEdgeFluxP3d{1}(:,:,1), InnerEdgeFluxS3d{1}(:,:,1) );
[ InnerSurface_frhs2d ] = edge2d.matEvaluateStrongFormEdgeRHS( InnerEdgeFluxM2d{1}(:,:,1), InnerEdgeFluxP2d{1}(:,:,1), InnerEdgeFluxS2d{1}(:,:,1) );

edge = mesh3d.BoundaryEdge;
edge2d = obj.mesh2d(1).BoundaryEdge;
[ BoundaryEdgefm3d{1}, BoundaryEdgefp3d{1} ] = edge.matEvaluateSurfValue( fphys3d );
[ BoundaryEdgefm3d{1}, BoundaryEdgefp3d{1} ] = obj.matImposeBoundaryCondition( edge, edge.nx, edge.ny, BoundaryEdgefm3d{1}, BoundaryEdgefp3d{1}, obj.fext3d{1} );
[ BoundaryFluxM3d{1} ] = obj.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, BoundaryEdgefm3d{1} );
[ BoundaryFluxM2d{1}(:,:,1) ] = edge.VerticalColumnIntegralField( BoundaryFluxM3d{1}(:,:,1) );
[ BoundaryEdgeFluxS3d{1} ] = obj.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, BoundaryEdgefm3d{1}(:,:,[4, 1, 2]), BoundaryEdgefp3d{1}(:,:,[4, 1, 2]), edge );
[ BoundaryEdgeFluxS2d{1} ] = edge.VerticalColumnIntegralField(BoundaryEdgeFluxS3d{1}(:,:,1) );

[ BoundarySurface_frhs3d ] = edge.matEvaluateStrongFormEdgeRHS( BoundaryFluxM3d{1}(:,:,1), BoundaryEdgeFluxS3d{1}(:,:,1) );
[ BoundarySurface_frhs2d ] = edge2d.matEvaluateStrongFormEdgeRHS( BoundaryFluxM2d{1}(:,:,1), BoundaryEdgeFluxS2d{1}(:,:,1) );









%     edge = mesh3d.InnerEdge;
%     edge2d = obj.mesh2d(1).InnerEdge;
%
%
%     InnerSurface_frhs3d = edge.matEvaluateStrongFromEdgeRHS( obj.InnerEdgeFluxM3d{1}(:,:,1), obj.InnerEdgeFluxP3d{1}(:,:,1), obj.InnerEdgeFluxS3d{1}(:,:,1) );
%     InnerSurface_frhs2d = edge2d.matEvaluateStrongFromEdgeRHS( obj.InnerEdgeFluxM2d{1}(:,:,1), obj.InnerEdgeFluxP2d{1}(:,:,1), obj.InnerEdgeFluxS2d{1}(:,:,1) );
%
%     edge = mesh3d.BoundaryEdge;
%     edge2d = obj.mesh2d(1).BoundaryEdge;
%
%     BoundarySurface_frhs3d = edge.matEvaluateStrongFormEdgeRHS( obj.BoundaryEdgeFluxM3d{1}(:,:,1), obj.BoundaryEdgeFluxS3d{1}(:,:,1) );
%     BoundarySurface_frhs2d = edge2d.matEvaluateStrongFromEdgeRHS( obj.BoundaryEdgeFluxM2d{1}(:,:,1), obj.BoundaryEdgeFluxS2d{1}(:,:,1) );

dHUHV =...
    ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,1)) + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,1)) ) - ...
    InnerSurface_frhs3d - BoundarySurface_frhs3d  + ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,2)) + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,2)) ) ;
% Term2d = ;
dHU2DHV2D = obj.meshUnion(1).Extend2dField( obj.mesh2d.rx .* (obj.mesh2d.cell.Dr * fphys2d{1}(:,:,2)) + obj.mesh2d.sx .* (obj.mesh2d.cell.Ds * fphys2d{1}(:,:,2) ) - ...
    InnerSurface_frhs2d - BoundarySurface_frhs2d + obj.mesh2d.ry .* (obj.mesh2d.cell.Dr * fphys2d{1}(:,:,3) ) + obj.mesh2d.sy .* (obj.mesh2d.cell.Ds * fphys2d{1}(:,:,3)) );

% edge = mesh3d.BottomEdge;
% fphys3d{1}(edge.GFToN1) = ( fphys3d{1}(edge.GFToN1) + fphys3d{1}(edge.GFToN2) )./2;
% fphys3d{1}(edge.GFToN2) = fphys3d{1}(edge.GFToN1);
%
% fphys3d{1}(edge.GFToN1 + mesh3d.cell.Np * mesh3d.K) = (fphys3d{1}(edge.GFToN1 + mesh3d.cell.Np * mesh3d.K)+...
%     fphys3d{1}(edge.GFToN2 + mesh3d.cell.Np * mesh3d.K))./2;
% fphys3d{1}(edge.GFToN2 + mesh3d.cell.Np * mesh3d.K) = fphys3d{1}(edge.GFToN1 + mesh3d.cell.Np * mesh3d.K);

% dHUHV =...
%     ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,1)) + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,1)) )  ...
%       + ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,2)) + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,2)) ) ;
% % Term2d = ;
% dHU2DHV2D = obj.meshUnion(1).Extend2dField( obj.mesh2d.rx .* (obj.mesh2d.cell.Dr * fphys2d{1}(:,:,2)) + obj.mesh2d.sx .* (obj.mesh2d.cell.Ds * fphys2d{1}(:,:,2) )  ...
%      + obj.mesh2d.ry .* (obj.mesh2d.cell.Dr * fphys2d{1}(:,:,3) ) + obj.mesh2d.sy .* (obj.mesh2d.cell.Ds * fphys2d{1}(:,:,3)) );

Omega =  mesh3d.VerticalIntegralField( -1 * dHUHV ) + dHU2DHV2D .* (1+mesh3d.z);

W = Omega + ( 1 + mesh3d.z ) .* ( - ( dHUHV ) - mesh3d.tz .* (mesh3d.cell.Dt * Omega) ) + ...
    fphys3d{1}(:,:,1)./fphys3d{1}(:,:,4) .*  mesh3d.z .* ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,4))...
    + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,4)) ) +  ...
    fphys3d{1}(:,:,1)./fphys3d{1}(:,:,4) .* ( mesh3d.rx .* (mesh3d.cell.Dr * fphys3d{1}(:,:,7))...
    + mesh3d.sx .* (mesh3d.cell.Ds * fphys3d{1}(:,:,7)) )+...
    fphys3d{1}(:,:,2)./fphys3d{1}(:,:,4) .* mesh3d.z.* ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,4))...
    + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,4)) )+...
    fphys3d{1}(:,:,2)./fphys3d{1}(:,:,4) .* ( mesh3d.ry .* (mesh3d.cell.Dr * fphys3d{1}(:,:,7))...
    + mesh3d.sy .* (mesh3d.cell.Ds * fphys3d{1}(:,:,7)) );
end

