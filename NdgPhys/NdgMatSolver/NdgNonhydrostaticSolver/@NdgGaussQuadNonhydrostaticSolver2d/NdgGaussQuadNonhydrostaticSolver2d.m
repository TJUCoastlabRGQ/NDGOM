classdef NdgGaussQuadNonhydrostaticSolver2d < NdgGaussQuadNonhydrostaticSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver2d( physClass, meshUnion )
            obj = obj@NdgGaussQuadNonhydrostaticSolver( physClass, meshUnion );
        end
    end
    
end

