classdef NdgQuadFreeWeakFormCentralVisSolver2d < NdgVisCentralFluxSolver &...
        NdgVisAuxiWeakFormSolver
    %NDGQUADFREEWEAKFORMCENTRALVISSOLVER2D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
         function obj = NdgQuadFreeWeakFormCentralVisSolver2d( phys )
            obj = obj@NdgVisCentralFluxSolver();
            obj = obj@NdgVisAuxiWeakFormSolver(phys);
        end       
    end
    
end

