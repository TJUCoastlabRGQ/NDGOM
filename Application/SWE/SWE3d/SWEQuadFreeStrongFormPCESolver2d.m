classdef SWEQuadFreeStrongFormPCESolver2d
    %SWEQUADFREESTRONGFORMPCESOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
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
                
                %Function used to calculate the two dimentional PCE volume term
                physClass.frhs2d{m} = -( ...
                    mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d{m}(:, :, 2) ) + ...
                    mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d{m}(:, :, 2) ) + ...
                    mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d{m}(:, :, 3) ) + ...
                    mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d{m}(:, :, 3) ) );
                
                % Function used to calculate the two dimentional PCE inner surface term
                InnerEdge = mesh2d.InnerEdge;
                [ fm, fp ] = InnerEdge.matEvaluateSurfValue( fphys );
                %> $\lambda = abs( max(sqrt{(gH^_)},sqrt{(gH^+)}))$
                lambda = abs( max( max( sqrt( obj.gra .* fm(:, :, 1) ), sqrt( obj.gra .* fp(:, :, 1) ) )) );
                FluxM = fm(:, :, 2) .* InnerEdge.nx + fm(:, :, 3) .* InnerEdge.ny;
                FluxP = fp(:, :, 2) .* InnerEdge.nx + fp(:, :, 3) .* InnerEdge.ny;
                %> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
                FluxS = 0.5 * ( FluxM + FluxP - bsxfun( @times, lambda, ( fp(:, :, 1) - fm(:, :, 1)  ) ) );
                physClass.frhs2d{m} = physClass.frhs2d{m} + InnerEdge.matEvaluateStrongFromEdgeRHS( FluxM, FluxP, FluxS );
                
                % Function used to calculate the two dimentional PCE boundary surface integration term
                BoundaryEdge = mesh2d.BoundaryEdge;
                [ fm, fp ] = BoundaryEdge.matEvaluateSurfValue( fphys );
                
                % apply clamped boundary condition
                ind = ( BoundaryEdge.ftype == enumBoundaryCondition.Clamped );
                fp(:, ind, 1) = fext(:, ind, 1);
                
                % apply slip wall boundary condition
                ind = ( BoundaryEdge.ftype == enumBoundaryCondition.SlipWall );
                Hun =  fm( :, ind, 2 ) .* BoundaryEdge.nx(:, ind) + fm( :, ind, 3).* BoundaryEdge.ny(:, ind);
                Hvn = -fm( :, ind, 2 ) .* BoundaryEdge.ny(:, ind) + fm( :, ind, 3).* BoundaryEdge.nx(:, ind);
                
                fp(:, ind, 2) = - Hun .* BoundaryEdge.nx(:, ind) - Hvn .* BoundaryEdge.ny(:, ind);
                fp(:, ind, 3) = - Hun .* BoundaryEdge.ny(:, ind) + Hvn .* BoundaryEdge.nx(:, ind);
                
                %> $\lambda = abs( max(sqrt{(gH^_)},sqrt{(gH^+)}))$
                lambda = abs( max( max( sqrt( obj.gra .* fm(:, :, 1) ), sqrt( obj.gra .* fp(:, :, 1) ) ) ) );
                % lambda = zeros(size(lambda));
                
                FluxM = fm(:, :, 2) .* BoundaryEdge.nx + fm(:, :, 3) .* BoundaryEdge.ny;
                FluxP = fp(:, :, 2) .* BoundaryEdge.nx + fp(:, :, 3) .* BoundaryEdge.ny;
                %> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
                FluxS = 0.5 * ( FluxM + FluxP - bsxfun( @times, lambda, ( fp(:, :, 1) - fm(:, :, 1) )) );
                
                physClass.frhs2d{m} = physClass.frhs2d{m} + BoundaryEdge.matEvaluateStrongFromEdgeRHS( FluxM, FluxS );
                
            end
        end
    end
    
end

