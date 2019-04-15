classdef NdgQuadFreeStrongCentralVisSolver3d < NdgAbstractVisSolver
    
    methods
        function obj = NdgQuadFreeStrongCentralVisSolver3d( phys, varId, rhsId )
            obj = obj@NdgAbstractVisSolver( phys, varId, rhsId );
        end
        
        function matEvaluateRHS( obj, fphys )
            matEvaluateAuxiVar( obj, fphys );
            matEvaluateOriVarRHS( obj, fphys );
        end
    end
    
    methods ( Access = protected )
        function matEvaluateAuxiVar( obj, fphys )
            matEvaluateAuxiVarVolumeKernel( obj, fphys );
            matEvaluateAuxiVarSurfaceKernel( obj, fphys );
            % todo: consider boundary condition
        end
        
        function matEvaluateAuxiVarVolumeKernel( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.phys.mesh3d(m);
                dfdr = mesh.cell.Dt * fphys{m}(:, :, obj.varId(1));
                dfds = mesh.cell.Dt * fphys{m}(:, :, obj.varId(2));
                obj.pzx{m} = mesh.tz .* dfdr;
                obj.pzy{m} = mesh.tz .* dfds;
            end
        end
        
        function matEvaluateAuxiVarSurfaceKernel( obj, fphys3d )
            for m = 1:obj.Nmesh
                mesh = obj.phys.mesh3d(m);
                edge3d = mesh.BottomEdge;
                [ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
                FluxM_1(:, :, 1) = edge3d.nz .* fm(:, :, 1);
                FluxM_1(:, :, 2) = edge3d.nz .* fm(:, :, 2);
                FluxP_1(:, :, 1) = edge3d.nz .* fp(:, :, 1);
                FluxP_1(:, :, 2) = edge3d.nz .* fp(:, :, 2);
                obj.pzx{m} = obj.pzx{m} - ...
                    edge3d.matEvaluateStrongFormEdgeCentralRHS( FluxM_1(:,:,1), FluxP_1(:,:,1) );
                obj.pzy{m} = obj.pzy{m} - ...
                    edge3d.matEvaluateStrongFormEdgeCentralRHS( FluxM_1(:,:,2), FluxP_1(:,:,2) );
                
                edge3d = mesh3d.BottomBoundaryEdge;
                [ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
                FluxM(:, :, 1) = edge3d.nz .* fm(:, :, 1);
                FluxM(:, :, 2) = edge3d.nz .* fm(:, :, 2);
                %> $|(Hu)^+ = (Hu)^-|_{\Omega = -1}$
                %> $|(Hv)^+ = (Hv)^-|_{\Omega = -1}$
                FluxS(:, :, 1) = edge3d.nz .* fp(:, :, 1);
                FluxS(:, :, 2) = edge3d.nz .* fp(:, :, 2);
                obj.pzx{m} = obj.pzx{m} ...
                    - edge3d.matEvaluateStrongFormEdgeRHS( FluxM(:,:,1), FluxS(:,:,1) );
                obj.pzy{m} = obj.pzy{m} ...
                    - edge3d.matEvaluateStrongFormEdgeRHS( FluxM(:,:,2), FluxS(:,:,2) );
                
                edge3d = mesh3d.SurfaceBoundaryEdge;
                [ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
                FluxM(:, :, 1) = edge3d.nz .* fm(:, :, 1);
                FluxM(:, :, 2) = edge3d.nz .* fm(:, :, 2);
                %> $|(Hu)^+ = (Hu)^-|_{\Omega = 0}$
                %> $|(Hv)^+ = (Hv)^-|_{\Omega = 0}$
                FluxS(:, :, 1) = edge3d.nz .* fp(:, :, 1);
                FluxS(:, :, 2) = edge3d.nz .* fp(:, :, 2);
                obj.pzx{m} = obj.pzx{m} ...
                    - edge3d.matEvaluateStrongFormEdgeRHS( FluxM(:,:,1), FluxS(:,:,1) );
                obj.pzy{m} = obj.pzy{m} ...
                    - edge3d.matEvaluateStrongFormEdgeRHS( FluxM(:,:,2), FluxS(:,:,2) );
            end
        end
        
        function matEvaluateOriVarRHS( obj, fphys )
            matEvaluateOriVarSurfaceKernel( obj, fphys );
            matEvaluateOriVarVolumeKernel( obj, fphys );
        end
        
        function matEvaluateOriVarVolumeKernel( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.phys.mesh3d(m);
                
                obj.phys.frhs{m}(:, :, 1) = obj.phys.frhs{m}(:, :, 1) + ...
                    mesh.rz .* (mesh.cell.Dr * obj.pzx{m}) + mesh.sz .* (mesh.cell.Ds * obj.pzx{m})...
                    + mesh.tz .* (mesh.cell.Dt * obj.pzx{m});
                obj.phys.frhs{m}(:, :, 2) = obj.phys.frhs{m}(:, :, 2) + ...
                    mesh.rz .* (mesh.cell.Dr * obj.pzy{m}) + mesh.sz .* (mesh.cell.Ds * obj.pzy{m})...
                    + mesh.tz .* (mesh.cell.Dt * obj.pzy{m});
            end
        end% func
        
        function matEvaluateOriVarSurfaceKernel( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.phys.mesh3d(m);
                obj.pzx{m} = obj.phys.miu{m} .*  obj.pzx{m};
                obj.pzy{m} = obj.phys.miu{m} .*  obj.pzy{m};
                
                edge = mesh.BottomEdge;
                [ fmx, fpx ] = edge.matEvaluateSurfValue( obj.pzx{m} );
                [ fmy, fpy ] = edge.matEvaluateSurfValue( obj.pzy{m} );
                fluxMx = edge.nz .* fmx; fluxMy = edge.nz .* fmy;
                fluxPx = edge.nz .* fpx; fluxPy = edge.nz .* fpy;
                obj.phys.frhs{m}(:, :, obj.rhsId{1}) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId{1}) + ...
                    edge.matEvaluateStrongFormEdgeCentralRHS( fluxMx, fluxPx );
                obj.phys.frhs{m}(:, :, obj.rhsId{2}) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId{2}) + ...
                    edge.matEvaluateStrongFormEdgeCentralRHS( fluxMy, fluxPy );
                
                % The wind term and the bottom friction term is considered is the source term fraction
                edge = mesh3d.SurfaceBoundaryEdge;
                [ fmx, ~ ] = edge.matEvaluateSurfValue( obj.pzx{m} );
                [ fmy, ~ ] = edge.matEvaluateSurfValue( obj.pzy{m} );
                fluxMx = edge.nz .* fmx; fluxMy = edge.nz .* fmy;
                fluxSx = zeros(size(fmx)); fluxSy = zeros(size(fmy));
                obj.phys.frhs{m}(:, :, obj.rhsId{1}) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId{1}) + ...
                    edge.matEvaluateStrongFormEdgeRHS( fluxMx, fluxSx );
                obj.phys.frhs{m}(:, :, obj.rhsId{2}) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId{2}) + ...
                    edge.matEvaluateStrongFormEdgeRHS( fluxMy, fluxSy );
                
                edge = mesh3d.BottomBoundaryEdge;
                [ fmx, ~ ] = edge.matEvaluateSurfValue( obj.pzx{m} );
                [ fmy, ~ ] = edge.matEvaluateSurfValue( obj.pzy{m} );
                fluxMx = edge.nz .* fmx; fluxMy = edge.nz .* fmy;
                fluxSx = zeros(size(fmx)); fluxSy = zeros(size(fmy));
                obj.phys.frhs{m}(:, :, obj.rhsId{1}) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId{1}) + ...
                    edge.matEvaluateStrongFormEdgeRHS( fluxMx, fluxSx );
                obj.phys.frhs{m}(:, :, obj.rhsId{2}) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId{2}) + ...
                    edge.matEvaluateStrongFormEdgeRHS( fluxMy, fluxSy );
                
            end
        end
    end
    
end

