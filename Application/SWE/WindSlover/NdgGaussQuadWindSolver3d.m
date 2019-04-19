classdef NdgGaussQuadWindSolver3d < AbstractWindTermSolver & ...
        NdgGaussQuadWeakFormSolver3d
    %NDGGAUSSQUADWINDSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadWindSolver3d( physObj, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver3d(physObj, meshUnion);
        end
        
        function evaluateWindTermRHS( obj, physClass, fphys )
            tempRHS = cell(physClass.Nmesh);
            for m = 1:physClass.Nmesh
                mesh3d = physClass.meshUnion(m);
                edge = mesh3d.SurfaceBoundaryEdge;
                %                [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                %                 fluxM = zeros(size(fm(:,:,1:2)));
                [ Taux, ~ ] = obj.matInterpolateToFaceGaussQuadraturePoint( obj.SBVfq, physClass.Taux{m}, physClass.Taux{m} );
                [ Tauy, ~ ] = obj.matInterpolateToFaceGaussQuadraturePoint( obj.SBVfq, physClass.Tauy{m}, physClass.Tauy{m} );
                
                Taux = obj.SBnz .* Taux;
                Tauy = obj.SBnz .* Tauy;
                
                tempRHS{m} = zeros(size(physObj.frhs{m}));
                
                tempRHS{m}(:,:,1) = obj.matAssembleSourceTermIntoRHS( edge, Taux, tempRHS{m}(:,:,1));
                tempRHS{m}(:,:,2) = obj.matAssembleSourceTermIntoRHS( edge, Tauy, tempRHS{m}(:,:,2));
                
                tempRHS{m} = permute( sum( bsxfun(@times, obj.invM{m}, ...
                    permute( permute( phys.frhs{m}, [1,3,2] ),[2,1,3] ) ), 2 ), [1,3,2]);
                
                physClass.frhs{m} = physClass.frhs{m} + tempRHS{m};
            end
        end
    end
    
end

