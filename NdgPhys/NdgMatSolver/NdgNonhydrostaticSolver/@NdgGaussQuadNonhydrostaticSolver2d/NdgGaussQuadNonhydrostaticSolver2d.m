classdef NdgGaussQuadNonhydrostaticSolver2d < NdgGaussQuadNonhydrostaticSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver2d( physClass, meshUnion )
            obj = obj@NdgGaussQuadNonhydrostaticSolver( physClass, meshUnion );
        end
    end
    
end

