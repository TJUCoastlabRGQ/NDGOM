classdef NdgQuadFreeStrongFormAdvSolver3d < NdgQuadFreeStrongFormSolver & ...
        NdgAbstractAdvSolver
    
    methods
        function obj = NdgQuadFreeStrongFormAdvSolver3d( phys )
            obj = obj@NdgQuadFreeStrongFormSolver( phys );
            obj = obj@NdgAbstractAdvSolver(phys);
        end
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, physClass, fphys )            
            % evaluate inner edge
            for m = 1:physClass.Nmesh
                mesh3d = physClass.meshUnion(m);
%                 mesh2d = physClass.mesh2d(m);                
                
                edge = mesh3d.InnerEdge;
%                 edge2d = mesh2d.InnerEdge;
                [ physClass.InnerEdgefm3d{m}, physClass.InnerEdgefp3d{m} ] = edge.matEvaluateSurfValue( fphys );
              
                [ physClass.InnerEdgeFluxM3d{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefm3d{m} );
%                 [ fluxM2d ] = edge.VerticalColumnIntegralField( fluxM(:,:,1) );
                [ physClass.InnerEdgeFluxP3d{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefp3d{m} );
%                 [ fluxP2d ] = edge.VerticalColumnIntegralField( fluxP(:,:,1) );
                [ physClass.InnerEdgeFluxS3d{m} ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, physClass.InnerEdgefm3d{m}(:,:,[4, 1, 2]), physClass.InnerEdgefp3d{m}(:,:,[4, 1, 2]), edge );
%                 [ fluxS2dOld ] = phys.matEvaluateSurfNumFlux( mesh2d, edge2d.nx, edge2d.ny, fm2d, fp2d, edge2d );
%                 [ fluxS2d ] = edge.VerticalColumnIntegralField( fluxS(:,:,1) );
                [ physClass.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( physClass.InnerEdgeFluxM3d{m}(:,:,[2,3]), physClass.InnerEdgeFluxP3d{m}(:,:,[2,3]), physClass.InnerEdgeFluxS3d{m}(:,:,[2,3]) );
%                 [ phys.frhs2d{m} ] = edge2d.matEvaluateStrongFromEdgeRHS( fluxM2d, fluxP2d, fluxS2d );

                edge = mesh3d.BoundaryEdge;
%                 edge2d = mesh2d.BoundaryEdge;
                [ physClass.BoundaryEdgefm3d{m}, physClass.BoundaryEdgefp3d{m} ] = edge.matEvaluateSurfValue( fphys );
    
                [ physClass.BoundaryEdgefm3d{m}, physClass.BoundaryEdgefp3d{m} ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, physClass.BoundaryEdgefm3d{m}, physClass.BoundaryEdgefp3d{m}, physClass.fext3d{m} );
                [ physClass.BoundaryEdgeFluxM3d{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm3d{m} );
%                 [ fluxM2d ] = edge.VerticalColumnIntegralField( fluxM(:,:,1) );
                [ physClass.BoundaryEdgeFluxS3d{m} ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, physClass.BoundaryEdgefm3d{m}(:,:,[4, 1, 2]), physClass.BoundaryEdgefp3d{m}(:,:,[4, 1, 2]), edge );
%                 [ fluxS2d ] = edge.VerticalColumnIntegralField( fluxS(:,:,1) );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( physClass.BoundaryEdgeFluxM3d{m}(:,:,[2,3]), physClass.BoundaryEdgeFluxS3d{m}(:,:,[2,3]) );
%                 [ phys.frhs2d{m} ] = phys.frhs2d{m} + edge2d.matEvaluateStrongFromEdgeRHS( fluxM2d, fluxS2d );
                
                % we note that for the three dimensional nonlinear shallow
                % water equation, Newmann boundary about the velocity is
                % imposed at both the surface and bottom boundary, under
                % this kind of circumastances, the numerical flux is the
                % same with the local numerical flux, so the contribution
                % of these term to the right hand is canceled for the
                % advective part
                
%                 edge = mesh3d.SurfaceBoundaryEdge;
%                 [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
%                 [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
%                 [ fluxS ] =  zeros(size(fluxM));
%                 [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
%                 
                edge = mesh3d.BottomEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );
%                 [ OmegafluxS(:,:,1) ] = edge.nz .* fm(:,:,1).*fm(:,:,3)./fm(:,:,4) .* ( edge.nz .* fm(:,:,3)>=0 ) + edge.nz .* fp(:,:,1).*fp(:,:,3)./fp(:,:,4) .* ( edge.nz .* fm(:,:,3)<0 );
%                 [ OmegafluxS(:,:,2) ] = edge.nz .* fm(:,:,2).*fm(:,:,3)./fm(:,:,4) .* ( edge.nz .* fm(:,:,3)>=0 ) + edge.nz .* fp(:,:,2).*fp(:,:,3)./fp(:,:,4) .* ( edge.nz .* fm(:,:,3)<0 );                
                
                [ OmegafluxS(:,:,1) ] = 0.5*edge.nz.*(fm(:,:,1).*fm(:,:,3)./fm(:,:,4)+fp(:,:,1).*fp(:,:,3)./fp(:,:,4));
                [ OmegafluxS(:,:,2) ] = 0.5*edge.nz.*(fm(:,:,2).*fm(:,:,3)./fm(:,:,4)+fp(:,:,2).*fp(:,:,3)./fp(:,:,4));
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), fluxP(:,:,[2,3]), OmegafluxS );
                
%                 edge = mesh3d.BottomBoundaryEdge;
%                 [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
%                 [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );            
% %                 [ fluxS ] =  zeros(size(fluxM));
%                 [ fluxS ] = fluxM;
%                 [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );    
            end
            
            for m = 1:physClass.Nmesh % calculate RHS term on each mesh
                mesh = physClass.meshUnion(m);
                [ E, G, H ] = physClass.matEvaluateFlux( mesh, fphys{m} );
                
                for i = 1:physClass.Nvar
%                     phys.frhs{m}(:,:,i) = ...
%                         phys.frhs{m}(:,:,i) + ...
%                         - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) ) ...
%                         - obj.sx{m}.*( obj.Ds{m} * E(:,:,i) ) ...
%                         - obj.tx{m}.*( obj.Dt{m} * E(:,:,i) ) ...
%                         - obj.ry{m}.*( obj.Dr{m} * G(:,:,i) ) ...
%                         - obj.sy{m}.*( obj.Ds{m} * G(:,:,i) ) ...
%                         - obj.ty{m}.*( obj.Dt{m} * G(:,:,i) ) ...
%                         - obj.rz{m}.*( obj.Dr{m} * H(:,:,i) ) ...
%                         - obj.sz{m}.*( obj.Ds{m} * H(:,:,i) ) ...
%                         - obj.tz{m}.*( obj.Dt{m} * H(:,:,i) );
                    physClass.frhs{m}(:,:,i) = ...
                        physClass.frhs{m}(:,:,i) ...
                        - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) ) ...
                        - obj.sx{m}.*( obj.Ds{m} * E(:,:,i) ) ...
                        - obj.ry{m}.*( obj.Dr{m} * G(:,:,i) ) ...
                        - obj.sy{m}.*( obj.Ds{m} * G(:,:,i) ) ...
                        - obj.tz{m}.*( obj.Dt{m} * H(:,:,i) );
                end
%                     phys.frhs2d{m}(:,:,1) = ...
%                         phys.frhs2d{m}(:,:,1) ...
%                         - mesh2d(m).rx.*( mesh2d(m).cell.Dr * fphys2d{m}(:,:,2) ) ...
%                         - mesh2d(m).sx.*( mesh2d(m).cell.Ds * fphys2d{m}(:,:,2)  ) ...
%                         - mesh2d(m).ry.*( mesh2d(m).cell.Dr * fphys2d{m}(:,:,3)  ) ...
%                         - mesh2d(m).sy.*( mesh2d(m).cell.Ds * fphys2d{m}(:,:,3)  );                
            end
        end
    end
    
end