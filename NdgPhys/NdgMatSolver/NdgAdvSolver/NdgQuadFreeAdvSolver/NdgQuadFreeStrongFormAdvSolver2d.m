classdef NdgQuadFreeStrongFormAdvSolver2d < NdgQuadFreeStrongFormSolver & ...
        NdgAbstractAdvSolver
    
    methods
        function obj = NdgQuadFreeStrongFormAdvSolver2d( phys )
            obj = obj@NdgQuadFreeStrongFormSolver( phys );
            obj = obj@NdgAbstractAdvSolver(phys);
        end
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, physClass, fphys )
            phys = obj.phys;
            
            % evaluate inner edge
            for m = 1:phys.Nmesh
                mesh = phys.meshUnion(m);                
                edge = mesh.InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
              
                [ fluxM ] = physClass.matEvaluateSurfFlux( mesh, edge.nx, edge.ny, fm, edge );
                [ fluxP ] = physClass.matEvaluateSurfFlux( mesh, edge.nx, edge.ny, fp, edge );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh, edge.nx, edge.ny, fm, fp, edge );
                [ physClass.frhs{m} ] = edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxP, fluxS );

                edge = mesh.BoundaryEdge;
                [ fm, physClass.BoundaryEdgefp{m} ] = edge.matEvaluateSurfValue( fphys );
    
                [ fm, physClass.BoundaryEdgefp{m} ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny, fm, physClass.BoundaryEdgefp{m}, physClass.fext{m} );
                [ fluxM ] = physClass.matEvaluateSurfFlux( mesh, edge.nx, edge.ny, fm, edge );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh, edge.nx, edge.ny, fm, physClass.BoundaryEdgefp{m}, edge );
                [ physClass.frhs{m} ] = physClass.frhs{m} + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );                
            end
            
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E, G ] = phys.matEvaluateFlux( mesh, fphys{m} );
%                 [ fm, fp ] = phys.matEvaluateSurfaceValue( mesh, fphys{m}, phys.fext{m} );
%                 [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.nx{m}, obj.ny{m}, fm, fp );
%                 [ flux ] = phys.matEvaluateSurfFlux( mesh, obj.nx{m}, obj.ny{m}, fm );
                
                for i = 1:phys.Nvar
                    phys.frhs{m}(:,:,i) = ...
                        phys.frhs{m}(:,:,i) + ...
                        - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) ) ...
                        - obj.sx{m}.*( obj.Ds{m} * E(:,:,i) ) ...
                        - obj.ry{m}.*( obj.Dr{m} * G(:,:,i) ) ...
                        - obj.sy{m}.*( obj.Ds{m} * G(:,:,i) ); ...
%                         + ( obj.LIFT{m} * ( obj.Js{m} .* ( flux(:,:,i) - fluxS(:,:,i) ) ))./ obj.J{m};
                end
            end
        end
    end
    
end

