classdef SWEGaussQuadWeakFormPCESolver2d < NdgGaussQuadWeakFormSolver
    %SWEQUADFREESTRONGFORMPCESOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        
        function  obj = SWEGaussQuadWeakFormPCESolver2d( phys, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( phys, meshUnion );
        end
        function evaluateAdvectionRHS( obj, physClass, fphys2d )
            
            % evaluate inner edge
                mesh2d = physClass.mesh2d(1);
                
                %  Function used to calculate the vertically averaged horizontal momentum term
                fq = obj.matInterpolateToVolumeGaussQuadraturePoint(obj.Vq{1}, fphys2d{1} );
                
                % Volume Integral
                [ physClass.frhs2d{1} ] = ...
                    + obj.Dr{1} * ( obj.rxwJ{1}.* (fq(:,:,2)) + obj.rywJ{1}.* ( fq(:,:,3) ) ) ...
                    + obj.Ds{1} * ( obj.sxwJ{1}.* (fq(:,:,2)) + obj.sywJ{1}.* ( fq(:,:,3) ) );
                
                % Function used to calculate the two dimentional PCE inner surface term
                InnerEdge = mesh2d.InnerEdge;
                [ fm, fp ] = InnerEdge.matEvaluateSurfValue( fphys2d );
                
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( InnerEdge, obj.IEFVfq{1}, fm, fp);
                
                %> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
                FluxS = physClass.matEvaluateSurfNumFlux( mesh2d, obj.IEnx{1}, obj.IEny{1}, fm, fp, InnerEdge );
                EdgeRHS = - ( obj.IELIFT{1} * ( obj.IEwJs{1} .* ( FluxS(:,:,1) ) ));
                physClass.frhs2d{1} = obj.matAssembleIntoRHS( InnerEdge, EdgeRHS, physClass.frhs2d{1});
                
                % Function used to calculate the two dimentional PCE boundary surface integration term
                BoundaryEdge = mesh2d.BoundaryEdge;
                [ fm, fp ] = BoundaryEdge.matEvaluateSurfValue( fphys2d );
                
                % apply clamped boundary condition
                ind = ( BoundaryEdge.ftype == enumBoundaryCondition.Clamped );
                fp(:, ind, 1) = physClass.fext2d{1}(:, ind, 1);
                fp(:, ind, 2) = physClass.fext2d{1}(:, ind, 2);
                fp(:, ind, 3) = physClass.fext2d{1}(:, ind, 3);
                % apply slip wall boundary condition
                ind = ( BoundaryEdge.ftype == enumBoundaryCondition.SlipWall );
                Hun =  fm( :, ind, 2 ) .* BoundaryEdge.nx(:, ind) + fm( :, ind, 3).* BoundaryEdge.ny(:, ind);
                Hvn = -fm( :, ind, 2 ) .* BoundaryEdge.ny(:, ind) + fm( :, ind, 3).* BoundaryEdge.nx(:, ind);
                
                fp(:, ind, 2) = - Hun .* BoundaryEdge.nx(:, ind) - Hvn .* BoundaryEdge.ny(:, ind);
                fp(:, ind, 3) = - Hun .* BoundaryEdge.ny(:, ind) + Hvn .* BoundaryEdge.nx(:, ind);
                
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( BoundaryEdge, obj.BEFVfq{1}, fm, fp);
                
                FluxS = physClass.matEvaluateSurfNumFlux( mesh2d, obj.BEnx{1}, obj.BEny{1}, fm, fp, BoundaryEdge );
                EdgeRHS = - ( obj.BELIFT{1} * ( obj.BEwJs{1} .* ( FluxS(:,:,1) ) ));
                
                physClass.frhs2d{1} = obj.matAssembleBoundaryAndSourceTermIntoRHS( BoundaryEdge, EdgeRHS, physClass.frhs2d{1});
                
                physClass.frhs2d{1} = permute( sum( ...
                    bsxfun(@times, obj.invM{1}, ...
                    permute( permute( physClass.frhs2d{1}, [1,3,2] ), ...
                    [2,1,3] ) ), 2 ), [1,3,2]);
        end
    end
    
end

