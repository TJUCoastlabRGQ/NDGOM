classdef SWEGaussQuadWeakFormPCESolver2d < NdgGaussQuadWeakFormSolver
    %SWEQUADFREESTRONGFORMPCESOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        
        function  obj = SWEGaussQuadWeakFormPCESolver2d( phys )
            obj = obj@NdgGaussQuadWeakFormSolver( phys );
        end
        function evaluateAdvectionRHS( obj, physClass, fphys2d, fphys3d, fext )
            
            fphys = cell(1);
            % evaluate inner edge
            for m = 1:physClass.Nmesh
                
                %  Function used to calculate the vertically averaged horizontal momentum term
                mesh3d = physClass.mesh3d(m);
                mesh2d = physClass.mesh2d(m);
                fphys2d{m}(:, :, 2) = mesh3d.VerticalColumnIntegralField( fphys3d{m}(:, :, 1) );
                fphys2d{m}(:, :, 3) = mesh3d.VerticalColumnIntegralField( fphys3d{m}(:, :, 2) );
                fphys{1} = fphys2d{m};
                fq = obj.matInterpolateToVolumeGaussQuadraturePoint(obj.Vq, fphys2d{m} );
                
                % Volume Integral
                 [ physClass.frhs2d{m} ] = ...
                        + obj.Dr{m} * ( obj.rxwJ{m}.* (fq(:,:,2)) + obj.rywJ{m}.* ( fq(:,:,3) ) ) ...
                        + obj.Ds{m} * ( obj.sxwJ{m}.* (fq(:,:,2)) + obj.sywJ{m}.* ( fq(:,:,3) ) );
                
                % Function used to calculate the two dimentional PCE inner surface term
                InnerEdge = mesh2d.InnerEdge;
                [ fm, fp ] = InnerEdge.matEvaluateSurfValue( fphys );
                
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( InnerEdge, obj.IEFVfq, fm, fp);
                %> $\lambda = abs( max(sqrt{(gH^_)},sqrt{(gH^+)}))$
                lambda = abs( max( max( sqrt( physClass.gra .* fm(:, :, 1) ), sqrt( physClass.gra .* fp(:, :, 1) ) )) );
                
                FluxM = fm(:, :, 2) .* obj.IEnx + fm(:, :, 3) .* obj.IEny;
                FluxP = fp(:, :, 2) .* obj.IEnx + fp(:, :, 3) .* obj.IEny;
                %> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
                FluxS = 0.5 * ( FluxM + FluxP - bsxfun( @times, lambda, ( fp(:, :, 1) - fm(:, :, 1)  ) ) );
                EdgeRHS = - ( obj.IELIFT{m} * ( obj.IEwJs{m} .* ( FluxS ) ));
                physClass.frhs2d{m} = obj.matAssembleIntoRHS( InnerEdge, EdgeRHS, physClass.frhs2d{m});
%                 physClass.frhs2d{m} = physClass.frhs2d{m} + InnerEdge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS );
                
                % Function used to calculate the two dimentional PCE boundary surface integration term
                BoundaryEdge = mesh2d.BoundaryEdge;
                [ fm, fp ] = BoundaryEdge.matEvaluateSurfValue( fphys );
                
                % apply clamped boundary condition
                ind = ( BoundaryEdge.ftype == enumBoundaryCondition.Clamped );
                fp(:, ind, 1) = fext{m}(:, ind, 1);
                
                % apply slip wall boundary condition
                ind = ( BoundaryEdge.ftype == enumBoundaryCondition.SlipWall );
                Hun =  fm( :, ind, 2 ) .* BoundaryEdge.nx(:, ind) + fm( :, ind, 3).* BoundaryEdge.ny(:, ind);
                Hvn = -fm( :, ind, 2 ) .* BoundaryEdge.ny(:, ind) + fm( :, ind, 3).* BoundaryEdge.nx(:, ind);
                
                fp(:, ind, 2) = - Hun .* BoundaryEdge.nx(:, ind) - Hvn .* BoundaryEdge.ny(:, ind);
                fp(:, ind, 3) = - Hun .* BoundaryEdge.ny(:, ind) + Hvn .* BoundaryEdge.nx(:, ind);
                
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( InnerEdge, obj.BEFVfq, fm, fp);
                
                %> $\lambda = abs( max(sqrt{(gH^_)},sqrt{(gH^+)}))$
                lambda = abs( max( max( sqrt( physClass.gra .* fm(:, :, 1) ), sqrt( physClass.gra .* fp(:, :, 1) ) ) ) );
                % lambda = zeros(size(lambda));
                
                FluxM = fm(:, :, 2) .* obj.BEnx + fm(:, :, 3) .* obj.BEny;
                FluxP = fp(:, :, 2) .* obj.BEnx + fp(:, :, 3) .* obj.BEny;
                %> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
                FluxS = 0.5 * ( FluxM + FluxP - bsxfun( @times, lambda, ( fp(:, :, 1) - fm(:, :, 1) )) );
                
                EdgeRHS = - ( obj.BELIFT{m} * ( obj.BEwJs{m} .* ( FluxS ) ));
                
                physClass.frhs2d{m} = obj.matAssembleIntoRHS( BoundaryEdge, EdgeRHS, physClass.frhs2d{m});
                
            end
        end
    end
    
end

