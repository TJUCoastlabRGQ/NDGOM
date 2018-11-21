classdef NdgVisAuxiAbstractSolver < handle
    %NDGVISAUXIABSTRACTSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
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

