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
                
                edge = mesh3d.InnerEdge;
                [ physClass.InnerEdgefm3d{m}, physClass.InnerEdgefp3d{m} ] = edge.matEvaluateSurfValue( fphys );
              
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefm3d{m} );
                [ fluxP ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefp3d{m} );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefm3d{m}, physClass.InnerEdgefp3d{m}, edge );
                [ physClass.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );

                edge = mesh3d.BoundaryEdge;
                [ physClass.BoundaryEdgefm3d{m}, ~ ] = edge.matEvaluateSurfValue( fphys );
    
                [ physClass.BoundaryEdgefm3d{m}, ~ ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm3d{m}, physClass.BoundaryEdgefp3d{m}, physClass.fext3d{m} );
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm3d{m} );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm3d{m}, physClass.BoundaryEdgefp3d{m}, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
        
                edge = mesh3d.BottomBoundaryEdge;
                [ physClass.BottomBoundaryEdgefm3d{m}, ~ ] = edge.matEvaluateSurfValue( fphys );
                [ physClass.BottomBoundaryEdgefm3d{m}, ~ ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, physClass.BottomBoundaryEdgefm3d{m}, physClass.BottomBoundaryEdgefp3d{m}, physClass.fext3d{m} );                
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.BottomBoundaryEdgefm3d{m} );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, physClass.BottomBoundaryEdgefm3d{m}, physClass.BottomBoundaryEdgefp3d{m}, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                edge = mesh3d.SurfaceBoundaryEdge;
                [ physClass.SurfaceBoundaryEdgefm3d{m}, ~ ] = edge.matEvaluateSurfValue( fphys );
                [ physClass.SurfaceBoundaryEdgefm3d{m}, ~ ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, physClass.SurfaceBoundaryEdgefm3d{m}, physClass.SurfaceBoundaryEdgefp3d{m}, physClass.fext3d{m} );                                
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.SurfaceBoundaryEdgefm3d{m} );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, physClass.SurfaceBoundaryEdgefm3d{m}, physClass.SurfaceBoundaryEdgefp3d{m}, edge );                
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                  
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