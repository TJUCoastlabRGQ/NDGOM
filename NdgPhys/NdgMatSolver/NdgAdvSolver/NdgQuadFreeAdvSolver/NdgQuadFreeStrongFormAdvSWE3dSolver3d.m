classdef NdgQuadFreeStrongFormAdvSWE3dSolver3d < NdgQuadFreeStrongFormAdvSolver3d
    
    methods
        function obj = NdgQuadFreeStrongFormAdvSWE3dSolver3d( phys )
            obj = obj@NdgQuadFreeStrongFormAdvSolver3d( phys );
        end
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, physClass, fphys )            
            % evaluate inner edge
            for m = 1:physClass.Nmesh
                mesh3d = physClass.meshUnion(m);
                
                edge = mesh3d.InnerEdge;
                [ physClass.InnerEdgefm{m}, physClass.InnerEdgefp{m} ] = edge.matEvaluateSurfValue( fphys );
              
                [ physClass.InnerEdgeFluxM{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefm{m} );
                [ physClass.InnerEdgeFluxP{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefp{m} );
                [ physClass.InnerEdgeFluxS{m} ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, physClass.InnerEdgefm{m}(:,:,[4, 1, 2]), physClass.InnerEdgefp{m}(:,:,[4, 1, 2]), edge );
                [ physClass.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( physClass.InnerEdgeFluxM{m}(:,:,[2,3]), physClass.InnerEdgeFluxP{m}(:,:,[2,3]), physClass.InnerEdgeFluxS{m}(:,:,[2,3]) );

                edge = mesh3d.BoundaryEdge;
                [ physClass.BoundaryEdgefm{m}, physClass.BoundaryEdgefp{m} ] = edge.matEvaluateSurfValue( fphys );
    
                [ physClass.BoundaryEdgefm{m}, physClass.BoundaryEdgefp{m} ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, physClass.BoundaryEdgefm{m}, physClass.BoundaryEdgefp{m}, physClass.fext3d{m} );
                [ physClass.BoundaryEdgeFluxM{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm{m} );
                [ physClass.BoundaryEdgeFluxS{m} ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, physClass.BoundaryEdgefm{m}(:,:,[4, 1, 2]), physClass.BoundaryEdgefp{m}(:,:,[4, 1, 2]), edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( physClass.BoundaryEdgeFluxM{m}(:,:,[2,3]), physClass.BoundaryEdgeFluxS{m}(:,:,[2,3]) );
                
        
                edge = mesh3d.BottomEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );                
                
                [ OmegafluxS(:,:,1) ] = edge.nz.* fp(:,:,1).*fp(:,:,3)./fp(:,:,4);
                [ OmegafluxS(:,:,2) ] = edge.nz.* fp(:,:,2).*fp(:,:,3)./fp(:,:,4);
%                 [ OmegafluxS(:,:,1) ] = 0.5*edge.nz.*( fm(:,:,1).*fm(:,:,3)./fm(:,:,4) + fp(:,:,1).*fp(:,:,3)./fp(:,:,4) );
%                 [ OmegafluxS(:,:,2) ] = 0.5*edge.nz.*( fm(:,:,2).*fm(:,:,3)./fm(:,:,4) + fp(:,:,2).*fp(:,:,3)./fp(:,:,4) );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), fluxP(:,:,[2,3]), OmegafluxS );
                
                edge = mesh3d.BottomBoundaryEdge;
                [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), zeros(size(fluxM(:,:,[2,3]) )));                

%                 edge = mesh3d.SurfaceBoundaryEdge;
%                 [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
%                 [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
%                 
%                 [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), physClass.SurfaceDate );                        
            end
            
            for m = 1:physClass.Nmesh % calculate RHS term on each mesh
                mesh = physClass.meshUnion(m);
                [ E, G, H ] = physClass.matEvaluateFlux( mesh, fphys{m} );
                
                for i = 1:physClass.Nvar
                    physClass.frhs{m}(:,:,i) = ...
                        physClass.frhs{m}(:,:,i) ...
                        - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) ) ...
                        - obj.sx{m}.*( obj.Ds{m} * E(:,:,i) ) ...
                        - obj.ry{m}.*( obj.Dr{m} * G(:,:,i) ) ...
                        - obj.sy{m}.*( obj.Ds{m} * G(:,:,i) ) ...
                        - obj.tz{m}.*( obj.Dt{m} * H(:,:,i) );
                end
             
            end
        end
    end
    
end