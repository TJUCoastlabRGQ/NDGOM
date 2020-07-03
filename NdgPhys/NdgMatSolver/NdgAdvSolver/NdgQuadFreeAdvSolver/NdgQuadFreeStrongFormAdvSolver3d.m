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
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
              
                [ fluxM ] = physClass.matEvaluateSurfFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, edge );
                [ fluxP ] = physClass.matEvaluateSurfFlux( mesh3d, edge.nx, edge.ny, edge.nz, fp, edge );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, fp, edge );
                [ physClass.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );

                edge = mesh3d.BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
    
                [ fm, fp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, fm, fp, physClass.fext{m} );
                [ fluxM ] = physClass.matEvaluateSurfFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, edge );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, fp, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                edge = mesh3d.BottomEdge;
                [ fm, fp ]  = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = physClass.matEvaluateSurfFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, edge );
                [ fluxP ] = physClass.matEvaluateSurfFlux( mesh3d, edge.nx, edge.ny, edge.nz, fp, edge );                
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, fp, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS );                

                edge = mesh3d.BottomBoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, fm, fp, physClass.BotEdgefext{m} );                
                [ fluxM ] = physClass.matEvaluateSurfFlux( edge, edge.nx, edge.ny, edge.nz, fm );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, fp, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                edge = mesh3d.SurfaceBoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, fm, fp, physClass.SurfEdgefext{m} );                                
                [ fluxM ] = physClass.matEvaluateSurfFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, edge );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, edge.nx, edge.ny, edge.nz, fm, fp, edge );                
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