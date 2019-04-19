classdef NdgGaussQuadFrictionSolver3d < AbstractFrictionTermSolver & ...
        NdgGaussQuadWeakFormSolver3d
    %NDGGAUSSQUADFRICTIONSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        
        function obj = NdgGaussQuadFrictionSolver3d( physObj, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver3d(physObj, meshUnion);
        end
        
        function evaluateFrictionTermRHS( obj, physObj, fphys )
            
            tempRHS = cell(physObj.Nmesh);
            for m = 1:physObj.Nmesh
                mesh3d = physObj.meshUnion(m);
                edge = mesh3d.BottomBoundaryEdge;
               [ fm, ~ ] = edge.matEvaluateSurfValue( fphys );
                fluxM = zeros(size(fm(:,:,1:2)));
               
                u = fm(:,:,1)./fm(:,:,4); v = fm(:,:,2)./fm(:,:,4); 
                Velocity = sqrt( u.^2 + v.^2 );
                friction(:,:,1) = physObj.Cf{m} .* u .* Velocity;
                friction(:,:,2) = physObj.Cf{m} .* v .* Velocity; 
                
                [ friction, ~ ] = obj.matInterpolateToFaceGaussQuadraturePoint( obj.BBFVfq, friction, friction );
                tempRHS{m} = zeros(size(physObj.frhs{m}));
                
                tempRHS{m}(:,:,1) = obj.matAssembleSourceTermIntoRHS( edge, friction(:,:,1), tempRHS{m}(:,:,1));
                tempRHS{m}(:,:,2) = obj.matAssembleSourceTermIntoRHS( edge, friction(:,:,2), tempRHS{m}(:,:,2));
                
                tempRHS{m} = permute( sum( bsxfun(@times, obj.invM{m}, ...
                    permute( permute( phys.frhs{m}, [1,3,2] ),[2,1,3] ) ), 2 ), [1,3,2]);  
                
              physObj.frhs{m} = physObj.frhs{m} + tempRHS{m};
                
                % frhs = frhs + cd*rouair*w10*windx/rouwater
%                 physObj.frhs{m}(:,:,2) = physObj.frhs{m}(:,:,2)...
%                     + edge.nz .*obj.Cf{m} .* v .* Velocity; 
                
            end
        end% func           
    end
    
end

