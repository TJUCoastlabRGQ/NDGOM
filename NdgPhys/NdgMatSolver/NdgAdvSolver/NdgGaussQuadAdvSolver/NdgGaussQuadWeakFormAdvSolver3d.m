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
                GQfphys = obj.matInterpolateToVolumeGaussQuadraturePoint( obj.Vq, fphys{m});
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
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.IEFVfq, fm, fp);
                
%                 [ fluxM ] = phys.matEvaluateSurfFlux( edge, obj.IEnx, obj.IEny, obj.IEnz, fm );
%                 [ fluxP ] = phys.matEvaluateSurfFlux( edge, obj.IEnx, obj.IEny, obj.IEnz, fp );
                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, obj.IEnx, obj.IEny, obj.IEnz, fm, fp );
                for i = 1:phys.Nvar
                    EdgeRHS = - ( obj.IELIFT{m} * ( obj.IEwJs{m} .* ( fluxS(:,:,i) ) ));
                    phys.frhs{m}(:,:,i) = obj.matAssembleIntoRHS( BoundaryEdge, EdgeRHS, phys.frhs{m}(:,:,i));
                end
                                
                edge = mesh3d.BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                % fext is not interpolated, wrong, need to be corrected
                [ fm, fp ] = phys.matImposeBoundaryCondition( edge, obj.BEnx, obj.BEny, obj.BEnz, fm, fp, phys.fext3d{m} );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.BEFVfq, fm, fp);

                [ fluxS ] = phys.matEvaluateSurfNumFlux( edge, edge.nx, edge.ny, edge.nz, fm, fp );
                for i = 1:phys.Nvar
                    EdgeRHS = - ( obj.IELIFT{m} * ( obj.IEwJs{m} .* ( fluxS(:,:,i) ) ));
                    phys.frhs{m}(:,:,i) = obj.matAssembleIntoRHS( BoundaryEdge, EdgeRHS, phys.frhs{m}(:,:,i));
                end
               
                phys.frhs{m} = permute( sum( ...
                    bsxfun(@times, obj.invM{m}, ...
                    permute( permute( phys.frhs{m}, [1,3,2] ), ...
                    [2,1,3] ) ), 2 ), [1,3,2]);                
            end
        end
    end
    
end

