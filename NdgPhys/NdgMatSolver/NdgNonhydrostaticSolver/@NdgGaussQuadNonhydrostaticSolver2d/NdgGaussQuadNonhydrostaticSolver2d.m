classdef NdgGaussQuadNonhydrostaticSolver2d < NdgNonhydrostaticSolver2d & ...
        NdgGaussQuadWeakFormSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver2d( physClass, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( physClass, meshUnion );
        end
    end
    
end

