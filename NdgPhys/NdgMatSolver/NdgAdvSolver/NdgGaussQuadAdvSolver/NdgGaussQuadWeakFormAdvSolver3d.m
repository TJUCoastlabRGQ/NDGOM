classdef NdgGaussQuadWeakFormAdvSolver3d < NdgGaussQuadWeakFormSolver3d & ...
        NdgAbstractAdvSolver
    %NDGGAUSSQUADWEAKFORMADVSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadWeakFormAdvSolver3d( phys, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver3d( phys, meshUnion );
            obj = obj@NdgAbstractAdvSolver(phys);
        end
        
        %> Call the flux subroutine from the NdgPhys object.
        function evaluateAdvectionRHS( obj, fphys )
            phys = obj.phys;
            
            for m = 1:phys.Nmesh % calculate RHS term on each mesh
                mesh = phys.meshUnion(m);
                GQfphys = obj.matInterpolateToVolumeGaussQuadraturePoint( obj.Vq{m}, fphys{m});
                [ E, G, H ] = phys.matEvaluateFlux( mesh, GQfphys );
                for i = 1:phys.Nvar
                    phys.frhs{m}(:,:,i) = ...
                        + obj.Dr{m} * ( obj.rxwJ{m} .* E(:,:,i) + obj.rywJ{m} .* G(:,:,i) + obj.rzwJ{m} .* H(:,:,i) ) ...
                        + obj.Ds{m} * ( obj.sxwJ{m} .* E(:,:,i) + obj.sywJ{m} .* G(:,:,i) + obj.szwJ{m} .* H(:,:,i) ) ...
                        + obj.Dt{m} * ( obj.txwJ{m} .* E(:,:,i) + obj.tywJ{m} .* G(:,:,i) + obj.tzwJ{m} .* H(:,:,i) );
                end
            end
            % evaluate inner edge
            for m = 1:phys.Nmesh
                mesh3d = phys.mesh3d(m);
                edge = mesh3d.InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.IEFVfq{m}, fm, fp);
                
                %                 [ fluxM ] = phys.matEvaluateSurfFlux( edge, obj.IEnx, obj.IEny, obj.IEnz, fm );
                %                 [ fluxP ] = phys.matEvaluateSurfFlux( edge, obj.IEnx, obj.IEny, obj.IEnz, fp );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, obj.IEnx{m}, obj.IEny{m}, obj.IEnz{m}, fm, fp );
                for i = 1:phys.Nvar
                    EdgeRHS = - ( obj.IELIFT{m} * ( obj.IEwJs{m} .* ( fluxS(:,:,i) ) ));
                    phys.frhs{m}(:,:,i) = obj.matAssembleIntoRHS( edge, EdgeRHS, phys.frhs{m}(:,:,i));
                end
                
                edge = mesh3d.BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                % fext is not interpolated, wrong, need to be corrected
                [ fm, fp ] = phys.matImposeBoundaryCondition( edge, edge.nx, edge.ny, edge.nz, fm, fp, phys.fext3d{m} );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.BEFVfq{m}, fm, fp);
                
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, obj.BEnx{m}, obj.BEny{m}, obj.BEnz{m}, fm, fp );
                for i = 1:phys.Nvar
                    EdgeRHS = - ( obj.IELIFT{m} * ( obj.BEwJs{m} .* ( fluxS(:,:,i) ) ));
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

