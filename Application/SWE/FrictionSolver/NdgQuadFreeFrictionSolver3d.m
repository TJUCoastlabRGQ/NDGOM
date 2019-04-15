classdef NdgQuadFreeFrictionSolver3d < AbstractFrictionTermSolver
    %NDGQUADFREEFRICTIONSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function evaluateFrictionTermRHS( obj, physObj, fphys )

            for m = 1:physObj.Nmesh
                mesh3d = physObj.meshUnion(m);
                edge = mesh3d.BottomBoundaryEdge;
               [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                u = fm(:,:,1)./fm(:,:,3); v = fm(:,:,2)./fm(:,:,3); 
                Velocity = sqrt( u.^2 + v.^2 );
%                 fluxS(:,:,1) = edge.nz .*obj.Cf{1} .* u .* Velocity;
%                 fluxS(:,:,2) = edge.nz .*obj.Cf{1} .* v .* Velocity;                
                
              physObj.frhs{m}(:,:,1) = physObj.frhs{m}(:,:,1)...
                    + edge.nz .*obj.Cf{m} .* u .* Velocity;
                
                % frhs = frhs + cd*rouair*w10*windx/rouwater
                physObj.frhs{m}(:,:,2) = physObj.frhs{m}(:,:,2)...
                    + edge.nz .*obj.Cf{m} .* v .* Velocity; 
                
            end
        end% func        
    end
    
end

