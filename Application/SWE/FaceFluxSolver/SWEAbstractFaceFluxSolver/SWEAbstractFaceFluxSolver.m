classdef SWEAbstractFaceFluxSolver
    %SWEABSTRACTFACEFLUXSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods(Abstract)
        [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, fm)
    end
    
end

