classdef SWEQuadFreeStrongFormPCESolver2d
    %SWEQUADFREESTRONGFORMPCESOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function evaluateAdvectionRHS( obj, physClass, fphys2d )
             % evaluate inner edge
            for m = 1:physClass.Nmesh
                
                %  Function used to calculate the vertically averaged horizontal momentum term
                mesh3d = physClass.meshUnion(m);
                mesh2d = physClass.mesh2d(m);
%                 fphys2d{m}(:, :, 2) = mesh3d.VerticalColumnIntegralField( fphys3d{m}(:, :, 1) );
%                 fphys2d{m}(:, :, 3) = mesh3d.VerticalColumnIntegralField( fphys3d{m}(:, :, 2) );                
                %Function used to calculate the two dimentional PCE volume term
                physClass.frhs2d{m} = -( ...
                    mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d{m}(:, :, 2) ) + ...
                    mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d{m}(:, :, 2) ) + ...
                    mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d{m}(:, :, 3) ) + ...
                    mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d{m}(:, :, 3) ) );
                
                % Function used to calculate the two dimentional PCE inner surface term
                InnerEdge = mesh3d.InnerEdge;
                InnerEdge2d = mesh2d.InnerEdge;
                physClass.InnerEdgeFluxM2d{m} = InnerEdge.VerticalColumnIntegralField( physClass.InnerEdgeFluxM3d{m}(:,:,1) );
                physClass.InnerEdgeFluxP2d{m} = InnerEdge.VerticalColumnIntegralField( physClass.InnerEdgeFluxP3d{m}(:,:,1) );
                %> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
                physClass.InnerEdgeFluxS2d{m} = InnerEdge.VerticalColumnIntegralField( physClass.InnerEdgeFluxS3d{m}(:,:,1) );
                physClass.frhs2d{m} = physClass.frhs2d{m} + InnerEdge2d.matEvaluateStrongFromEdgeRHS( physClass.InnerEdgeFluxM2d{m}, physClass.InnerEdgeFluxP2d{m}, physClass.InnerEdgeFluxS2d{m} );
                
                % Function used to calculate the two dimentional PCE boundary surface integration term
                BoundaryEdge = mesh3d.BoundaryEdge;
                BoundaryEdge2d = mesh2d.BoundaryEdge;

%                 FluxP = fp(:, :, 2) .* BoundaryEdge.nx + fp(:, :, 3) .* BoundaryEdge.ny;
                %> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(H^+ - H^-)$
%                 FluxS = 0.5 * ( FluxM + FluxP - bsxfun( @times, lambda, ( fp(:, :, 1) - fm(:, :, 1) )) );
                physClass.BoundaryEdgeFluxM2d{m} = BoundaryEdge.VerticalColumnIntegralField( physClass.BoundaryEdgeFluxM3d{m}(:,:,1) );
                physClass.BoundaryEdgeFluxS2d{m} = BoundaryEdge.VerticalColumnIntegralField( physClass.BoundaryEdgeFluxS3d{m}(:,:,1) );
                
                physClass.frhs2d{m} = physClass.frhs2d{m} + BoundaryEdge2d.matEvaluateStrongFromEdgeRHS( physClass.BoundaryEdgeFluxM2d{m}, physClass.BoundaryEdgeFluxS2d{m} );
                
            end
        end
    end
    
end

