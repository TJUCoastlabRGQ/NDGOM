classdef NdghydrostaticSolver2d
    %NDGHYDROSTATICSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdghydrostaticSolver2d(PhysClass)
            %             obj = obj@NdgAbstractNonhydrostaticSolver(PhysClass);
        end
        
        function fphys = NdgNonhydrostaticUpdata(obj, PhysClass, fphys)
            % doing nothing
        end
        
        function fphys = NdgConservativeNonhydrostaticUpdata(obj,PhysClass, fphys, rk, intRK, dt)
            % doing nothing
        end
        
    end
    
    
end

