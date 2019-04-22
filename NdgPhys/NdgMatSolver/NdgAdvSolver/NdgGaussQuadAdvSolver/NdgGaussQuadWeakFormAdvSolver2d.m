classdef NdgGaussQuadWeakFormAdvSolver2d < NdgGaussQuadWeakFormSolver & NdgAbstractAdvSolver
    
    methods
        function obj = NdgGaussQuadWeakFormAdvSolver2d( phys, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( phys, meshUnion );
            obj = obj@NdgAbstractAdvSolver( phys );
        end
        
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                fq = zeros( mesh.cell.Nq, mesh.K, phys.Nfield );
                for i = 1:phys.Nfield
                    fq(:,:,i) = obj.Vq{m} * fphys{m}(:,:,i);
                end
                [ E, G ] = phys.matEvaluateFlux( mesh, fq );
                
                for i = 1:phys.Nvar
                    phys.frhs{m}(:,:,i) = ...
                        + obj.Dr{m} * ( obj.rxwJ{m} .* E(:,:,i) + obj.rywJ{m} .* G(:,:,i)) ...
                        + obj.Ds{m} * ( obj.sxwJ{m} .* E(:,:,i) + obj.sywJ{m} .* G(:,:,i)) ...
                        + obj.Dt{m} * ( obj.txwJ{m} .* E(:,:,i) + obj.tywJ{m} .* G(:,:,i));
                end
                
                edge = mesh.InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.IEFVfq{m}, fm, fp);
                
                [ fluxS ] = phys.matEvaluateSurfNumFlux( mesh, obj.IEnx{m}, obj.IEny{m}, fm, fp );
                
                for i = 1:phys.Nvar
                    EdgeRHS = - ( obj.IELIFT{m} * ( obj.IEwJs{m} .* ( fluxS(:,:,i) ) ));
                    phys.frhs{m}(:,:,i) = obj.matAssembleIntoRHS( edge, EdgeRHS, phys.frhs{m}(:,:,i));
                end
                
                edge = mesh.BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                % fext is not interpolated, wrong, need to be corrected,
                % for the 2D version, the boundary condition is not considered at present
                [ fm, fp ] = phys.matImposeBoundaryCondition( edge, edge.nx, edge.ny, fm, fp, phys.fext );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.BEFVfq{m}, fm, fp);
                
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, obj.BEnx{m}, obj.BEny{m}, fm, fp );
                for i = 1:phys.Nvar
                    EdgeRHS = - ( obj.BELIFT{m} * ( obj.BEwJs{m} .* ( fluxS(:,:,i) ) ));
                    phys.frhs{m}(:,:,i) = obj.matAssembleBoundaryAndSourceTermIntoRHS( edge, EdgeRHS, phys.frhs{m}(:,:,i));
                end
                
                
                for i = 1:phys.Nvar
                    phys.frhs{m}(:,:,i) = permute( sum( ...
                        bsxfun(@times, obj.invM{m}, ...
                        permute( permute( phys.frhs{m}(:,:,i), [1,3,2] ), ...
                        [2,1,3] ) ), 2 ), [1,3,2]);
                end
            end
        end
    end
    
end

