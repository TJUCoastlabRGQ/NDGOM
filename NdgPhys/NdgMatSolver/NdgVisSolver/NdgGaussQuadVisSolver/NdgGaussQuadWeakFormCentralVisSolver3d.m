classdef NdgGaussQuadWeakFormCentralVisSolver3d < NdgGaussQuadWeakFormSolver3d & ...
        NdgAbstractVisSolver
    %NDGGAUSSQUADWEAKFORMVISSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        tempRHS
    end
    
    methods
        function obj = NdgGaussQuadWeakFormCentralVisSolver3d( phys, varId, rhsId )
            obj = obj@NdgAbstractVisSolver( phys, varId, rhsId );
            obj = obj@NdgGaussQuadWeakFormSolver3d( phys, phys.meshUnion );
            obj.tempRHS = cell(phys.Nmesh);
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
                fq = obj.matInterpolateToVolumeGaussQuadraturePoint(obj.Vq, fphys{m} );
                
                % Volume Integral
                obj.pzx{m} = obj.Dt{m} * ( obj.tzwJ{m}.* (fq(:,:,obj.varId(1))));
                obj.pzy{m} = obj.Dt{m} * ( obj.tzwJ{m}.* (fq(:,:,obj.varId(2))));
            end
        end
        
        function matEvaluateAuxiVarSurfaceKernel( obj, fphys3d )
            for m = 1:obj.Nmesh
                mesh = obj.phys.mesh3d(m);
                edge3d = mesh.BottomEdge;
                [ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge3d, obj.BOTFVfq, fm, fp);
                FluxM_1(:, :, 1) = obj.BOTnz .* fm(:, :, 1);
                FluxM_1(:, :, 2) = obj.BOTnz .* fm(:, :, 2);
                FluxP_1(:, :, 1) = obj.BOTnz .* fp(:, :, 1);
                FluxP_1(:, :, 2) = obj.BOTnz .* fp(:, :, 2);
                FluxS = ( FluxM_1 + FluxP_1 )/2;
                EdgeRHSPzx = ( obj.BOTLIFT{m} * ( obj.BOTwJs{m} .* ( FluxS(:,:,1) ) ));
                EdgeRHSPzy = ( obj.BOTLIFT{m} * ( obj.BOTwJs{m} .* ( FluxS(:,:,2) ) ));
                
                obj.pzx{m} = obj.matAssembleIntoRHS( edge3d, EdgeRHSPzx, obj.pzx{m});
                obj.pzy{m} = obj.matAssembleIntoRHS( edge3d, EdgeRHSPzy, obj.pzy{m});
                
                edge3d = mesh.BottomBoundaryEdge;
                [ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge3d, obj.BBFVfq, fm, fp);
                FluxM(:, :, 1) = obj.BBnz .* fm(:, :, 1);
                FluxM(:, :, 2) = obj.BBnz .* fm(:, :, 2);
                %> $|(Hu)^+ = (Hu)^-|_{\Omega = -1}$
                %> $|(Hv)^+ = (Hv)^-|_{\Omega = -1}$
                FluxS(:, :, 1) = obj.BBnz .* fp(:, :, 1);
                FluxS(:, :, 2) = obj.BBnz .* fp(:, :, 2);
                EdgeRHSPzx = ( obj.BBLIFT{m} * ( obj.BBwJs{m} .* ( FluxS(:,:,1) ) ));
                EdgeRHSPzy = ( obj.BBLIFT{m} * ( obj.BBwJs{m} .* ( FluxS(:,:,2) ) ));
                
                obj.pzx{m} = obj.matAssembleIntoRHS( edge3d, EdgeRHSPzx, obj.pzx{m});
                obj.pzy{m} = obj.matAssembleIntoRHS( edge3d, EdgeRHSPzy, obj.pzy{m});
                
                edge3d = mesh.SurfaceBoundaryEdge;
                [ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
                [ fm, fp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge3d, obj.BBFVfq, fm, fp);
                FluxM(:, :, 1) = obj.SBnz .* fm(:, :, 1);
                FluxM(:, :, 2) = obj.SBnz .* fm(:, :, 2);
                %> $|(Hu)^+ = (Hu)^-|_{\Omega = 0}$
                %> $|(Hv)^+ = (Hv)^-|_{\Omega = 0}$
                FluxS(:, :, 1) = obj.SBnz .* fp(:, :, 1);
                FluxS(:, :, 2) = obj.SBnz .* fp(:, :, 2);
                EdgeRHSPzx = ( obj.SBLIFT{m} * ( obj.SBwJs{m} .* ( FluxS(:,:,1) ) ));
                EdgeRHSPzy = ( obj.SBLIFT{m} * ( obj.SBwJs{m} .* ( FluxS(:,:,2) ) ));
                
                obj.pzx{m} = obj.matAssembleIntoRHS( edge3d, EdgeRHSPzx, obj.pzx{m});
                obj.pzy{m} = obj.matAssembleIntoRHS( edge3d, EdgeRHSPzy, obj.pzy{m});
                
                obj.pzx{m} = permute( sum( bsxfun(@times, obj.invM{m}, ...
                    permute( permute( obj.pzx{m}, [1,3,2] ),[2,1,3] ) ), 2 ), [1,3,2]);
                obj.pzy{m} = permute( sum(bsxfun(@times, obj.invM{m}, ...
                    permute( permute( obj.pzy{m}, [1,3,2] ),[2,1,3] ) ), 2 ), [1,3,2]);
            end
        end
        
        function matEvaluateOriVarRHS( obj, fphys )
            matEvaluateOriVarSurfaceKernel( obj, fphys );
            matEvaluateOriVarVolumeKernel( obj, fphys );
        end
        
        function matEvaluateOriVarVolumeKernel( obj, fphys )
            for m = 1:obj.Nmesh
%                 mesh = obj.phys.mesh3d(m);
                obj.pzx{m} = obj.phys.miu{m} .*  obj.pzx{m};
                obj.pzy{m} = obj.phys.miu{m} .*  obj.pzy{m};
                
                pzx = obj.matInterpolateToVolumeGaussQuadraturePoint(obj.Vq, obj.pzx{m} );
                pzy = obj.matInterpolateToVolumeGaussQuadraturePoint(obj.Vq, obj.pzy{m} );
                
                obj.tempRHS{m}(:,:,1) = obj.Dt{m} * ( obj.tzwJ{m}.* pzx);
                obj.tempRHS{m}(:,:,2) = obj.Dt{m} * ( obj.tzwJ{m}.* pzy);
            end
        end% func          
        
        function matEvaluateOriVarSurfaceKernel( obj, fphys )
            for m = 1:obj.Nmesh
                mesh = obj.phys.mesh3d(m);
                
                %only one mesh considered
                edge = mesh.BottomEdge;
                [ fmx, fpx ] = edge.matEvaluateSurfValue( obj.pzx );
                [ fmx, fpx ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge3d, obj.BOTFVfq, fmx, fpx );
                [ fmy, fpy ] = edge.matEvaluateSurfValue( obj.pzy );
                [ fmy, fpy ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge3d, obj.BOTFVfq, fmy, fpy );
                
                fluxMx = obj.BOTnz .* fmx; fluxMy = obj.BOTnz .* fmy;
                fluxPx = obj.BOTnz .* fpx; fluxPy = obj.BOTnz .* fpy;
                
                FluxS(:, :, 1) = (fluxMx + fluxPx)./2;
                FluxS(:, :, 2) = (fluxMy + fluxPy)./2;
                EdgeRHSPzx = ( obj.BOTLIFT{m} * ( obj.BOTwJs{m} .* ( FluxS(:,:,1) ) ));
                EdgeRHSPzy = ( obj.BOTLIFT{m} * ( obj.BOTwJs{m} .* ( FluxS(:,:,2) ) ));
                obj.tempRHS{m}(:,:,1) = obj.matAssembleIntoRHS( edge3d, EdgeRHSPzx, obj.tempRHS{m}(:,:,1));
                obj.tempRHS{m}(:,:,2) = obj.matAssembleIntoRHS( edge3d, EdgeRHSPzy, obj.tempRHS{m}(:,:,2));
                
                obj.phys.frhs{m}(:, :, obj.rhsId(1))  = obj.phys.frhs{m}(:, :, obj.rhsId(1)) + ...
                    permute( sum( bsxfun(@times, obj.invM{m}, ...
                    permute( permute( obj.tempRHS{m}(:,:,1), [1,3,2] ),[2,1,3] ) ), 2 ), [1,3,2]);

                obj.phys.frhs{m}(:, :, obj.rhsId(2))  = obj.phys.frhs{m}(:, :, obj.rhsId(2)) + ...
                    permute( sum( bsxfun(@times, obj.invM{m}, ...
                    permute( permute( obj.tempRHS{m}(:,:,2), [1,3,2] ),[2,1,3] ) ), 2 ), [1,3,2]);                
                
                % The wind term and the bottom friction term is considered is the source term fraction, So the bottom and surface boundary contribution to the RHS here is zero.   
            end
        end      
    end
end

