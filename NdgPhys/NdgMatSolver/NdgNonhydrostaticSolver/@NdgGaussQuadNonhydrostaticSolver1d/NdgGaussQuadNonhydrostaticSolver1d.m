classdef NdgGaussQuadNonhydrostaticSolver1d < NdgAbstractNonhydrostaticSolver & ...
        NdgGaussQuadWeakFormSolver
        
    %NDGGAUSSQUADNONHYDROSTATICSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgGaussQuadNonhydrostaticSolver1d( physClass, meshUnion )
            obj = obj@NdgGaussQuadWeakFormSolver( physClass, meshUnion );
        end
    end
    
end

