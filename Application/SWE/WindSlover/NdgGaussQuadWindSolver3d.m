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
            
            for m = 1:physClass.Nmesh
                mesh3d = physClass.meshUnion(m);
                edge = mesh3d.SurfaceBoundaryEdge;
               [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                fluxM = zeros(size(fm(:,:,1:2)));
               
               fluxS(:,:,1) = edge.nz .*physClass.Taux{m}; 
               fluxS(:,:,2) = edge.nz .*physClass.Tauy{m};
                
              physClass.frhs{m} = physClass.frhs{m}...
                    + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );   
            end
        end        
    end
    
end

