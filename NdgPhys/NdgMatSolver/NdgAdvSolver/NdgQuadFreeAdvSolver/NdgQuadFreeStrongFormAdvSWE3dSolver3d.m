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
                [ physClass.InnerEdgefm3d{m}, physClass.InnerEdgefp3d{m} ] = edge.matEvaluateSurfValue( fphys );
              
                [ physClass.InnerEdgeFluxM3d{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefm3d{m} );
                [ physClass.InnerEdgeFluxP3d{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefp3d{m} );
                [ physClass.InnerEdgeFluxS3d{m} ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, physClass.InnerEdgefm3d{m}(:,:,[4, 1, 2]), physClass.InnerEdgefp3d{m}(:,:,[4, 1, 2]), edge );
                [ physClass.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( physClass.InnerEdgeFluxM3d{m}(:,:,[2,3]), physClass.InnerEdgeFluxP3d{m}(:,:,[2,3]), physClass.InnerEdgeFluxS3d{m}(:,:,[2,3]) );

                edge = mesh3d.BoundaryEdge;
                [ physClass.BoundaryEdgefm3d{m}, physClass.BoundaryEdgefp3d{m} ] = edge.matEvaluateSurfValue( fphys );
    
                [ physClass.BoundaryEdgefm3d{m}, physClass.BoundaryEdgefp3d{m} ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, physClass.BoundaryEdgefm3d{m}, physClass.BoundaryEdgefp3d{m}, physClass.fext3d{m} );
                [ physClass.BoundaryEdgeFluxM3d{m} ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm3d{m} );
                [ physClass.BoundaryEdgeFluxS3d{m} ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, physClass.BoundaryEdgefm3d{m}(:,:,[4, 1, 2]), physClass.BoundaryEdgefp3d{m}(:,:,[4, 1, 2]), edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( physClass.BoundaryEdgeFluxM3d{m}(:,:,[2,3]), physClass.BoundaryEdgeFluxS3d{m}(:,:,[2,3]) );
                
        
                edge = mesh3d.BottomEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );                
                
                [ OmegafluxS(:,:,1) ] = 0.5*edge.nz.*(fm(:,:,1).*fm(:,:,3)./fm(:,:,4)+fp(:,:,1).*fp(:,:,3)./fp(:,:,4));
                [ OmegafluxS(:,:,2) ] = 0.5*edge.nz.*(fm(:,:,2).*fm(:,:,3)./fm(:,:,4)+fp(:,:,2).*fp(:,:,3)./fp(:,:,4));
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM(:,:,[2,3]), fluxP(:,:,[2,3]), OmegafluxS );
                
                
                  
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