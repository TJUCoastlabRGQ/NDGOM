classdef NdgGaussQuadNonhydrostaticSolver1d < NdgAbstractNonhydrostaticSolver & ...
        NdgGaussQuadWeakFormSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver1d( physClass, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( physClass, meshUnion );
        end
    end
    
end

