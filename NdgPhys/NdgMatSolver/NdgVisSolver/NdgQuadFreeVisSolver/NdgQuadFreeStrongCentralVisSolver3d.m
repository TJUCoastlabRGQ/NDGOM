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
                for fld = 1:obj.Nfield
                    id = obj.rhsId(fld);
                    
                    obj.phys.frhs{m}(:, :, id) = obj.phys.frhs{m}(:, :, id) + ...
                        mesh.rz .* (mesh.cell.Dr * fphys{m}(:,:,3+id)) + mesh.sz .* (mesh.cell.Ds * fphys{m}(:,:,3+id))...
                        + mesh.tz .* (mesh.cell.Dt * fphys{m}(:,:,3+id));
                end
            end
        end% func
        
        function matEvaluateOriVarSurfaceKernel( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.phys.mesh3d(m);
                fphys{m}(:,:,4) = obj.phys.miu{m} .*  obj.pzx{m};
                fphys{m}(:,:,5) = obj.phys.miu{m} .*  obj.pzy{m};
                
                edge = mesh.BottomEdge;
                [ fmx, fpx ] = edge.matEvaluateSurfValue( fphys{m}(:,:,4) );
                [ fmy, fpy ] = edge.matEvaluateSurfValue( fphys{m}(:,:,5) );
                fluxMx = edge.nz .* fmx; fluxMy = edge.nz .* fmy;
                fluxPx = edge.nz .* fpx; fluxPy = edge.nz .* fpy;
                obj.phys.frhs{m}(:, :, obj.rhsId{1}) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId{1}) + ...
                    edge.matEvaluateStrongFormEdgeCentralRHS( fluxMx, fluxPx );
                
                obj.phys.frhs{m}(:, :, obj.rhsId{2}) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId{2}) + ...
                    edge.matEvaluateStrongFormEdgeCentralRHS( fluxMy, fluxPy );
                % The wind term and the bottom friction term is considered is the source term fraction
            end
        end
    end
    
end

