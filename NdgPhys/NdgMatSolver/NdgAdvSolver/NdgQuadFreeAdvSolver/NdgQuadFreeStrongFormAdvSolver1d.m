classdef NdgQuadFreeStrongFormAdvSolver1d < NdgAbstractAdvSolver ...
        & NdgQuadFreeStrongFormSolver
    
    methods
        function obj = NdgQuadFreeStrongFormAdvSolver1d( phys )
            obj = obj@NdgQuadFreeStrongFormSolver( phys );
            obj = obj@NdgAbstractAdvSolver(phys);
        end
        
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                edge = phys.meshUnion(m).InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( mesh, edge.nx, fm, edge );
                [ fluxP ] = phys.matEvaluateSurfFlux( mesh, edge.nx, fp, edge );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, edge.nx, fm, fp, edge );
                [ phys.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );
                
                edge = phys.meshUnion(m).BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = phys.matImposeBoundaryCondition( edge, edge.nx, fm, fp, phys.fext{m} );
                [ fluxM ] = phys.matEvaluateSurfFlux( mesh, edge.nx, fm, edge );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, edge.nx, fm, fp, edge );
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxS );
                
            end
            
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                [ E ] = phys.matEvaluateFlux( mesh, fphys{m} );          
                for i = 1:phys.Nvar
                    phys.frhs{m}(:,:,i) = ...
                        phys.frhs{m}(:,:,i) + ...
                        - obj.rx{m}.*( obj.Dr{m} * E(:,:,i) );
                end
            end
        end
    end
    
end

