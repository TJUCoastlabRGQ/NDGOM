classdef NdgQuadFreeWindSolver3d < AbstractWindTermSolver
    %NDGQUADFREEWINDSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgQuadFreeWindSolver3d()
            %doing nothing
        end
        
        function evaluateWindTermRHS( obj, physClass, fphys )
            
%             fluxS(:,:,1) = edge.nz .*obj.Taux{1}; fluxS(:,:,2) = edge.nz .*obj.Tauy{1};
            for m = 1:physClass.Nmesh
                mesh3d = physClass.meshUnion(m);
                edge = mesh3d.SurfaceBoundaryEdge;
              physClass.frhs{m}(:,:,1) = physClass.frhs{m}(:,:,1)...
                    + edge.nz .*obj.Taux{1};
                
                % frhs = frhs + cd*rouair*w10*windx/rouwater
                physClass.frhs{m}(:,:,2) = physClass.frhs{m}(:,:,2)...
                    + edge.nz .*obj.Tauy{1};
                
            end
        end
    end
    
end

