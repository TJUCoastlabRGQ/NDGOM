classdef NdgQuadFreeStrongFormCentralVisSolver2d < NdgVisCentralFluxSolver &...
        NdgVisAuxiStrongFormSolver
    %NDGQUADFREESTRONGFORMCENTRALVISSOLVER2D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgQuadFreeStrongFormCentralVisSolver2d( phys )
            obj = obj@NdgVisCentralFluxSolver();
            obj = obj@NdgVisAuxiStrongFormSolver(phys);
        end
    end
    
end

