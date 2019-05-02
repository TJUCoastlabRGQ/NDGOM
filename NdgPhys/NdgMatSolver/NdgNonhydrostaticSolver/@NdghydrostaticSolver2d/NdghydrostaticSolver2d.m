classdef NdghydrostaticSolver2d
    %NDGHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
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

