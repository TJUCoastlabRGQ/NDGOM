classdef NdgQuadFreeStrongFormCentralVisSolver2d < NdgVisCentralFluxSolver &...
        NdgVisAuxiStrongFormSolver
    %NDGQUADFREESTRONGFORMCENTRALVISSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = NdgQuadFreeStrongFormCentralVisSolver2d( phys )
            obj = obj@NdgVisCentralFluxSolver();
            obj = obj@NdgVisAuxiStrongFormSolver(phys);
        end
    end
    
end

