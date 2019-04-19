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
                u = fm(:,:,1)./fm(:,:,4); v = fm(:,:,2)./fm(:,:,4);
                Velocity = sqrt( u.^2 + v.^2 );
                friction(:,:,1) = physObj.Cf{m} .* u .* Velocity;
                friction(:,:,2) = physObj.Cf{m} .* v .* Velocity;
                
                [ friction, ~ ] = obj.matInterpolateToFaceGaussQuadraturePoint(  edge, obj.BBFVfq{m}, friction, friction );
                
                % we note that this part is added as a supplement of the
                % diffusion operator, so the direction vector should be included
                friction(:,:,1) = obj.BBnz{m} .* friction(:,:,1);
                friction(:,:,2) = obj.BBnz{m} .* friction(:,:,2);
                
                TEMPFriction(:,:,1) = ( obj.BBLIFT{m} * ( obj.BBwJs{m} .* ( friction(:,:,1) ) ));
                TEMPFriction(:,:,2) = ( obj.BBLIFT{m} * ( obj.BBwJs{m} .* ( friction(:,:,2) ) ));
                
                tempRHS{m} = zeros(size(physObj.frhs{m}));
                %                 %                 ( obj, edge, FRHS, RHS)
                tempRHS{m}(:,:,1) = obj.matAssembleSourceTermIntoRHS( edge, TEMPFriction(:,:,1), tempRHS{m}(:,:,1));
                tempRHS{m}(:,:,2) = obj.matAssembleSourceTermIntoRHS( edge, TEMPFriction(:,:,2), tempRHS{m}(:,:,2));
                
                tempRHS{m}(:,:,1) = permute( sum( bsxfun(@times, obj.invM{m}, ...
                    permute( permute( tempRHS{m}(:,:,1), [1,3,2] ),[2,1,3] ) ), 2 ), [1,3,2]);
                tempRHS{m}(:,:,2) = permute( sum( bsxfun(@times, obj.invM{m}, ...
                    permute( permute( tempRHS{m}(:,:,2), [1,3,2] ),[2,1,3] ) ), 2 ), [1,3,2]);
                
                %                 physObj.frhs{m}(:,:,1) = physObj.frhs{m}(:,:,1) + obj.matAssembleSourceTermIntoRHS( edge, friction(:,:,1), physObj.frhs{m}(:,:,1));
                %                 physObj.frhs{m}(:,:,2) = physObj.frhs{m}(:,:,2)  + obj.matAssembleSourceTermIntoRHS( edge, friction(:,:,2),physObj.frhs{m}(:,:,2));
                % frhs = frhs + cd*rouair*w10*windx/rouwater
                physObj.frhs{m} = physObj.frhs{m} + tempRHS{m};
                
            end
        end% func
    end
    
end

