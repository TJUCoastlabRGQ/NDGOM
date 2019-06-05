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
                edge = phys.meshUnion(m).InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, fm );
                [ fluxP ] = phys.matEvaluateSurfFlux( edge, edge.nx, fp );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, edge.nx, fm, fp );
                [ phys.frhs{m} ] = edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxP, fluxS );
                
                edge = phys.meshUnion(m).BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
%                 [ fm, fp ] = phys.matImposeBoundaryCondition( edge, edge.nx, fm, fp, phys.fext{m} );
                [ fluxM ] = phys.matEvaluateSurfFlux( edge, edge.nx, fm );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, edge.nx, fm, fp );
                [ phys.frhs{m} ] = phys.frhs{m} + edge.matEvaluateStrongFromEdgeRHS( fluxM, fluxS );
                
                
                %                 mesh = phys.meshUnion(m);
                %                 [ E ] = phys.matEvaluateFlux( mesh, fphys{m} );
                %                 [ fm, fp ] = phys.matEvaluateSurfaceValue( mesh, fphys{m}, phys.fext{m} );
                %                 [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.nx{m}, fm, fp );
                %                 [ flux ] = phys.matEvaluateSurfFlux( mesh, obj.nx{m}, fm );
                %
                %                 for i = 1:phys.Nvar
                %                     [ phys.frhs{m}(:,:,i) ] = ...
                %                         - mesh.rx.*( mesh.cell.Dr * E(:,:,i) ) ...
                %                         + ( obj.LIFT{m} * ( obj.Js{m} .* ( flux(:,:,i) - fluxS(:,:,i) ) ))./ obj.J{m};
                %                 end
                
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

