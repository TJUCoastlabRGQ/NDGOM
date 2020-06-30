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
                [ physClass.InnerEdgefm{m}, physClass.InnerEdgefp{m} ] = edge.matEvaluateSurfValue( fphys );
              
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefm{m} );
                [ fluxP ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefp{m} );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, physClass.InnerEdgefm{m}, physClass.InnerEdgefp{m}, edge );
                [ physClass.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );

                edge = mesh3d.BoundaryEdge;
                [ physClass.BoundaryEdgefm{m}, ~ ] = edge.matEvaluateSurfValue( fphys );
    
                [ physClass.BoundaryEdgefm{m}, ~ ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm{m}, physClass.BoundaryEdgefp{m}, physClass.fext3d{m} );
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm{m} );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, physClass.BoundaryEdgefm{m}, physClass.BoundaryEdgefp{m}, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                edge = mesh3d.BottomEdge;
                [ fm, fp ]  = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxP ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fp );                
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, fp, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS );                

                edge = mesh3d.BottomBoundaryEdge;
                [ physClass.BottomBoundaryEdgefm{m}, ~ ] = edge.matEvaluateSurfValue( fphys );
                [ physClass.BottomBoundaryEdgefm{m}, ~ ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, physClass.BottomBoundaryEdgefm{m}, physClass.BottomBoundaryEdgefp{m}, physClass.fext3d{m} );                
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.BottomBoundaryEdgefm{m} );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, physClass.BottomBoundaryEdgefm{m}, physClass.BottomBoundaryEdgefp{m}, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                edge = mesh3d.SurfaceBoundaryEdge;
                [ physClass.SurfaceBoundaryEdgefm{m}, ~ ] = edge.matEvaluateSurfValue( fphys );
                [ physClass.SurfaceBoundaryEdgefm{m}, ~ ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, physClass.SurfaceBoundaryEdgefm{m}, physClass.SurfaceBoundaryEdgefp{m}, physClass.fext3d{m} );                                
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, physClass.SurfaceBoundaryEdgefm{m} );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, physClass.SurfaceBoundaryEdgefm{m}, physClass.SurfaceBoundaryEdgefp{m}, edge );                
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