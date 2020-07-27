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
        function evaluateAdvectionRHS( obj, physClass, fphys )
            
            for m = 1:physClass.Nmesh % calculate RHS term on each mesh
                mesh = physClass.meshUnion(m);
                GQfphys = obj.matInterpolateToVolumeGaussQuadraturePoint( obj.Vq{m}, fphys{m});
                [ E, G, H ] = physClass.matEvaluateFlux( mesh, GQfphys );
                for i = 1:physClass.Nvar
                    physClass.frhs{m}(:,:,i) = ...
                        + obj.Dr{m} * ( obj.rxwJ{m} .* E(:,:,i) + obj.rywJ{m} .* G(:,:,i) + obj.rzwJ{m} .* H(:,:,i) ) ...
                        + obj.Ds{m} * ( obj.sxwJ{m} .* E(:,:,i) + obj.sywJ{m} .* G(:,:,i) + obj.szwJ{m} .* H(:,:,i) ) ...
                        + obj.Dt{m} * ( obj.txwJ{m} .* E(:,:,i) + obj.tywJ{m} .* G(:,:,i) + obj.tzwJ{m} .* H(:,:,i) );
                end
            end
            % evaluate inner edge
            for m = 1:physClass.Nmesh
                mesh3d = physClass.meshUnion(m);
                edge = mesh3d.InnerEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.IEFVfq{m}, fm, fp);
                
                %                 [ fluxM ] = phys.matEvaluateSurfFlux( edge, obj.IEnx, obj.IEny, obj.IEnz, fm );
                %                 [ fluxP ] = phys.matEvaluateSurfFlux( edge, obj.IEnx, obj.IEny, obj.IEnz, fp );
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, obj.IEnx{m}, obj.IEny{m}, fm(:,:,[4 1 2]), fp(:,:,[4 1 2]), edge );
                for i = 1:physClass.Nvar
                    %because only the last two numerical flux returned by
                    %HLL Riemann solver is used here
                    EdgeRHS = - ( obj.IELIFT{m} * ( obj.IEwJs{m} .* ( fluxS(:,:,i + 1 ) ) ));
                    physClass.frhs{m}(:,:,i) = obj.matAssembleIntoRHS( edge, EdgeRHS, physClass.frhs{m}(:,:,i));
                end
                
                edge = mesh3d.BoundaryEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                % fext is not interpolated, wrong, need to be corrected
                [ fm, fp ] = physClass.matImposeBoundaryCondition( edge, edge.nx, edge.ny,  fm, fp, physClass.fext3d{m} );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.BEFVfq{m}, fm, fp);
                
                [ fluxS ] = physClass.matEvaluateSurfNumFlux( mesh3d, obj.BEnx{m}, obj.BEny{m}, fm(:,:,[4 1 2]), fp(:,:,[4 1 2]), edge );
                for i = 1:physClass.Nvar
                    EdgeRHS = - ( obj.BELIFT{m} * ( obj.BEwJs{m} .* ( fluxS(:,:,i+1) ) ));
                    physClass.frhs{m}(:,:,i) = obj.matAssembleBoundaryAndSourceTermIntoRHS( edge, EdgeRHS, physClass.frhs{m}(:,:,i));
                end

                edge = mesh3d.BottomEdge;
                [ fm, fp ] = edge.matEvaluateSurfValue( fphys );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.BOTFVfq{m}, fm, fp);
                [ OmegafluxS(:,:,1) ] = obj.BOTnz{m} .* fm(:,:,1).*fm(:,:,3)./fm(:,:,4) .* ( obj.BOTnz{m} .* fm(:,:,3)>=0 ) + obj.BOTnz{m} .* fp(:,:,1).*fp(:,:,3)./fp(:,:,4) .* ( obj.BOTnz{m} .* fm(:,:,3)<0 );
                [ OmegafluxS(:,:,2) ] = obj.BOTnz{m} .* fm(:,:,2).*fm(:,:,3)./fm(:,:,4) .* ( obj.BOTnz{m} .* fm(:,:,3)>=0 ) + obj.BOTnz{m} .* fp(:,:,2).*fp(:,:,3)./fp(:,:,4) .* ( obj.BOTnz{m} .* fm(:,:,3)<0 );               
                for i = 1:physClass.Nvar
                    EdgeRHS = - ( obj.BOTLIFT{m} * ( obj.BOTwJs{m} .* ( OmegafluxS(:,:,i) ) ));
                    physClass.frhs{m}(:,:,i) = obj.matAssembleIntoRHS( edge, EdgeRHS, physClass.frhs{m}(:,:,i));
                end         
                
                %> This part is added for validation purpose, generally, this part is not needed
                edge = mesh3d.SurfaceBoundaryEdge;
                [omegafluxS(:,:,1), omegafluxS(:,:,2)] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, obj.SBFVfq{m}, physClass.SurfaceDate(:,:,1), physClass.SurfaceDate(:,:,2));
                for i = 1:physClass.Nvar
                    EdgeRHS = - ( obj.SBLIFT{m} * ( obj.SBwJs{m} .* ( omegafluxS(:,:,i) ) ));
                    physClass.frhs{m}(:,:,i) = obj.matAssembleBoundaryAndSourceTermIntoRHS( edge, EdgeRHS, physClass.frhs{m}(:,:,i));
                end                
                
                for i = 1:physClass.Nvar
                    physClass.frhs{m}(:,:,i) = permute( sum( ...
                        bsxfun(@times, obj.invM{m}, ...
                        permute( permute( physClass.frhs{m}(:,:,i), [1,3,2] ), ...
                        [2,1,3] ) ), 2 ), [1,3,2]);
                end
            end
        end
    end
    
end

