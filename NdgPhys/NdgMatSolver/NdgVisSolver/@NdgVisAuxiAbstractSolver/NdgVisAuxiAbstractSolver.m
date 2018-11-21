classdef NdgVisAuxiAbstractSolver < handle
    %NDGVISAUXIABSTRACTSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties( SetAccess = protected )
        phys
    end
    
    methods
        function obj = NdgVisAuxiAbstractSolver( phys )
            obj.phys = phys;
        end
    end
    
    methods(Abstract)
        evaluateViscosityRHS( obj, fphys );
    end
    
end

