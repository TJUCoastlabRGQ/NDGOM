classdef NdgQuadFreeStrongCentralVisSolver3d < NdgAbstractVisSolver
    % At present, only the vertical diffusion term is considered, both the
    % wind term and the friction term need to be added here and not the
    % source term part
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
        end
        
        function matEvaluateAuxiVarVolumeKernel( obj, fphys )
          %> This part is used to calculated the volume contribution to terms '$ pzx = \frac{\partial Hu}{\partial \sigma} and pzy = \frac{\partial Hv}{\partial \sigma}$'  
            for m = 1:obj.Nmesh
                mesh = obj.phys.meshUnion(m);
                dfdr = mesh.cell.Dt * fphys{m}(:, :, obj.varId(1));
                dfds = mesh.cell.Dt * fphys{m}(:, :, obj.varId(2));
                obj.pzx{m} = mesh.tz .* dfdr;
                obj.pzy{m} = mesh.tz .* dfds;
            end
        end
        
        function matEvaluateAuxiVarSurfaceKernel( obj, fphys3d )
            %> This part is used to calculated the surface contribution to terms '$ pzx = \frac{\partial Hu}{\partial \sigma} and pzy = \frac{\partial Hv}{\partial \sigma}$'
            for m = 1:obj.Nmesh
                mesh = obj.phys.meshUnion(m);
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
                
%                 edge3d = mesh.BottomBoundaryEdge;
%                 [ fm, ~ ] = edge3d.matEvaluateSurfValue( fphys3d );
                
                % Actually, this is a Newmann boundary, for such a boundary
                % the exterior value is the same with the inner boundaries,
                % and the local flux term and the numerical flux term
                % canceled here. 
                
%                 FluxM(:, :, 1) = edge3d.nz .* fm(:, :, 1);
%                 FluxM(:, :, 2) = edge3d.nz .* fm(:, :, 2);
%                 %> $|(Hu)^+ = (Hu)^-|_{\Omega = -1}$
%                 %> $|(Hv)^+ = (Hv)^-|_{\Omega = -1}$
%                 u = fm(:,:,1)./fm(:,:,4); v = fm(:,:,2)./fm(:,:,4); 
%                 Velocity = sqrt( u.^2 + v.^2 );
                
                
%                 fluxS(:,:,1) = edge3d.nz .*obj.phys.Cf{m} .* u .* Velocity;
%                 fluxS(:,:,2) = edge3d.nz .*obj.phys.Cf{m} .* v .* Velocity; 
%                 
%               physObj.frhs{m} = physObj.frhs{m}...
%                     + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
%                 
%                 
%                 
%                 FluxS(:, :, 1) = edge3d.nz .* fp(:, :, 1);
%                 FluxS(:, :, 2) = edge3d.nz .* fp(:, :, 2);
%                 obj.pzx{m} = obj.pzx{m} ...
%                     - edge3d.matEvaluateStrongFormEdgeRHS( FluxM(:,:,1), FluxS(:,:,1) );
%                 obj.pzy{m} = obj.pzy{m} ...
%                     - edge3d.matEvaluateStrongFormEdgeRHS( FluxM(:,:,2), FluxS(:,:,2) );

                % Actually, this is a Newmann boundary, for such a boundary
                % the exterior value is the same with the inner boundaries,
                % and the local flux term and the numerical flux term
                % canceled here. 
                
%                 edge3d = mesh.SurfaceBoundaryEdge;
%                 [ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
%                 FluxM(:, :, 1) = edge3d.nz .* fm(:, :, 1);
%                 FluxM(:, :, 2) = edge3d.nz .* fm(:, :, 2);
%                 %> $|(Hu)^+ = (Hu)^-|_{\Omega = 0}$
%                 %> $|(Hv)^+ = (Hv)^-|_{\Omega = 0}$
%                 FluxS(:, :, 1) = edge3d.nz .* fp(:, :, 1);
%                 FluxS(:, :, 2) = edge3d.nz .* fp(:, :, 2);
%                 obj.pzx{m} = obj.pzx{m} ...
%                     - edge3d.matEvaluateStrongFormEdgeRHS( FluxM(:,:,1), FluxS(:,:,1) );
%                 obj.pzy{m} = obj.pzy{m} ...
%                     - edge3d.matEvaluateStrongFormEdgeRHS( FluxM(:,:,2), FluxS(:,:,2) );
            end
        end
        
        function matEvaluateOriVarRHS( obj, fphys )
            matEvaluateOriVarSurfaceKernel( obj, fphys );
            matEvaluateOriVarVolumeKernel( obj, fphys );
        end
        
        function matEvaluateOriVarVolumeKernel( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.phys.meshUnion(m);
                
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
                mesh = obj.phys.meshUnion(m);
                %> Before here, pzx is equal to '$\frac{\partial hu}{\partial \sigma}$' and pzy '$\frac{\partial hv}{\partial \sigma}$'
                %> for the subsequent calculation, they are firstly multiplied by '$\frac{\mu}{H^2}$'
                obj.pzx{m} = fphys{m}(:,:,5)./(fphys{m}(:,:,4)).^2 .*  obj.pzx{m};
                obj.pzy{m} = fphys{m}(:,:,5)./(fphys{m}(:,:,4)).^2 .*  obj.pzy{m};
                
                edge = mesh.BottomEdge;
                [ fmx, fpx ] = edge.matEvaluateSurfValue( obj.pzx );
                [ fmy, fpy ] = edge.matEvaluateSurfValue( obj.pzy );
                fluxMx = edge.nz .* fmx; fluxMy = edge.nz .* fmy;
                fluxPx = edge.nz .* fpx; fluxPy = edge.nz .* fpy;
                obj.phys.frhs{m}(:, :, obj.rhsId(1)) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId(1)) + ...
                    edge.matEvaluateStrongFormEdgeCentralRHS( fluxMx, fluxPx );
                obj.phys.frhs{m}(:, :, obj.rhsId(2)) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId(2)) + ...
                    edge.matEvaluateStrongFormEdgeCentralRHS( fluxMy, fluxPy );
                
                edge = mesh.SurfaceBoundaryEdge;
                [ fmx, ~ ] = edge.matEvaluateSurfValue( obj.pzx );
                [ fmy, ~ ] = edge.matEvaluateSurfValue( obj.pzy );
                fluxMx = edge.nz .* fmx; fluxMy = edge.nz .* fmy;
                %>Wind stress to be considered in this part. For this version, the wind stress at the 
                %> surface is imposed as '$\frac{\mu}{H^2}\frac{\partial Hu}{\partial sigma} = \tau_x$'
                %> and '$\frac{\mu}{H^2}\frac{\partial Hv}{\partial sigma} = \tau_y$'.
                %> Compared the wind stress in our study with that in FVCOM, the FVCOM version is recovered
                %> by multiplying the stess with \frac{1}{\rho_0}
                
                fluxSx = edge.nz .* obj.phys.WindTaux{m};
                fluxSy = edge.nz .* obj.phys.WindTauy{m};
                
                
                obj.phys.frhs{m}(:, :, obj.rhsId(1)) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId(1)) + ...
                    edge.matEvaluateStrongFormEdgeRHS( fluxMx, fluxSx );
                obj.phys.frhs{m}(:, :, obj.rhsId(2)) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId(2)) + ...
                    edge.matEvaluateStrongFormEdgeRHS( fluxMy, fluxSy );
                
                %only one mesh considered
                edge = mesh.BottomBoundaryEdge;
                [ fmx, ~ ] = edge.matEvaluateSurfValue( obj.pzx );
                [ fmy, ~ ] = edge.matEvaluateSurfValue( obj.pzy );
                fluxMx = edge.nz .* fmx; fluxMy = edge.nz .* fmy;
                [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                
                %> Friction to be considered in this part. For this version, the friction at the bottom is given as
                %> '$\frac{\partial u}{\partial sigma} = \frac{H}{\miu}C_f u\sqrt(u^2+v^2)$' and 
                %> '$\frac{\partial v}{\partial sigma} = \frac{H}{\miu}C_f v\sqrt(u^2+v^2)$'.
                %> it's noted that '$C_f = \frac{C_d}{\rho_0}$' when we set the current version
                %> to be the same with the FVCOM. With all of those condition, we know that
                 %> '$\frac{\partial Hu}{\partial sigma} = \frac{H^2}{\miu}C_f u\sqrt(u^2+v^2)$' and 
                %> '$\frac{\partial Hv}{\partial sigma} = \frac{H^2}{\miu}C_f v\sqrt(u^2+v^2)$'. We multiply this by \
                %> '$\frac{\mu}{H^2}$' to arrive at '$\frac{\mu}{H^2}\frac{\partial Hu}{\partial sigma} = C_f u\sqrt(u^2+v^2)$'
                %> and '$\frac{\mu}{H^2}\frac{\partial Hv}{\partial sigma} = C_f v\sqrt(u^2+v^2)$'
                u = fm(:,:,1)./fm(:,:,4); v = fm(:,:,2)./fm(:,:,4); 
                Velocity = sqrt( u.^2 + v.^2 );                
                fluxSx = edge.nz .* obj.phys.Cf{m}.* u .* Velocity;
                fluxSy = edge.nz .* obj.phys.Cf{m}.* v .* Velocity;

                obj.phys.frhs{m}(:, :, obj.rhsId(1)) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId(1)) + ...
                    edge.matEvaluateStrongFormEdgeRHS( fluxMx, fluxSx );
                obj.phys.frhs{m}(:, :, obj.rhsId(2)) = ...
                    obj.phys.frhs{m}(:, :, obj.rhsId(2)) + ...
                    edge.matEvaluateStrongFormEdgeRHS( fluxMy, fluxSy );
                              
                
            end
        end
    end
    
end

