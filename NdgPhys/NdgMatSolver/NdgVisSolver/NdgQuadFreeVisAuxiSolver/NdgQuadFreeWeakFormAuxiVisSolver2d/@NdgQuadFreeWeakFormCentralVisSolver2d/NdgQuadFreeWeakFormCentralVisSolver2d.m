classdef NdgQuadFreeWeakFormCentralVisSolver2d < NdgVisCentralFluxSolver &...
        NdgVisAuxiWeakFormSolver
    %NDGQUADFREEWEAKFORMCENTRALVISSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
         function obj = NdgQuadFreeWeakFormCentralVisSolver2d( phys )
            obj = obj@NdgVisCentralFluxSolver();
            obj = obj@NdgVisAuxiWeakFormSolver(phys);
        end       
    end
    
end

