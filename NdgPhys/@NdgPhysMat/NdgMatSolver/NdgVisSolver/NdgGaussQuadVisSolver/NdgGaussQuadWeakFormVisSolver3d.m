classdef NdgGaussQuadWeakFormVisSolver3d < NdgGaussQuadWeakFormSolver3d & ...
         NdgAbstractVisSolver
    %NDGGAUSSQUADWEAKFORMVISSOLVER3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadWeakFormVisSolver3d( phys, varId, rhsId )
            obj = obj@NdgAbstractVisSolver( phys, varId, rhsId );
            obj = obj@NdgGaussQuadWeakFormCentralVisSolver3d( phys, phys.meshUnion );
        end
        
%         function matEvaluateRHS( obj, fphys )
%             matEvaluateAuxiVar( obj, fphys );
%             matEvaluateOriVarRHS( obj, fphys );
%         end
    end
    
end

