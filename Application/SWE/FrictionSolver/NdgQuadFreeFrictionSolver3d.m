classdef NdgQuadFreeFrictionSolver3d < AbstractFrictionTermSolver
    %NDGQUADFREEFRICTIONSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function evaluateFrictionTermRHS( obj, physObj, fphys )
% This part need to be deleted, as friction is not a source term but a
% boundary condition here
            for m = 1:physObj.Nmesh
                mesh3d = physObj.meshUnion(m);
                edge = mesh3d.BottomBoundaryEdge;
               [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                fluxM = zeros(size(fm(:,:,1:2)));
               
                u = fm(:,:,1)./fm(:,:,4); v = fm(:,:,2)./fm(:,:,4); 
                Velocity = sqrt( u.^2 + v.^2 );
                fluxS(:,:,1) = edge.nz .*physObj.Cf{m} .* u .* Velocity;
                fluxS(:,:,2) = edge.nz .*physObj.Cf{m} .* v .* Velocity; 
                
              physObj.frhs{m} = physObj.frhs{m}...
                    + edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
                
                % frhs = frhs + cd*rouair*w10*windx/rouwater
%                 physObj.frhs{m}(:,:,2) = physObj.frhs{m}(:,:,2)...
%                     + edge.nz .*obj.Cf{m} .* v .* Velocity; 
                
            end
        end% func        
    end
    
end

