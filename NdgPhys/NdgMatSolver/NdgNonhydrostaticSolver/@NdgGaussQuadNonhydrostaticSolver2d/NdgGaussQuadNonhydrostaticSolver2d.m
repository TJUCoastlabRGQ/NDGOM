classdef NdgGaussQuadNonhydrostaticSolver2d < NdgNonhydrostaticSolver2d & ...
        NdgGaussQuadWeakFormSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver2d( physClass, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( physClass, meshUnion );
        end
    end
    
end

