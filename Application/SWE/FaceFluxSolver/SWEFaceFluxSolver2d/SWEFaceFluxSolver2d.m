classdef SWEFaceFluxSolver2d < SWEAbstractFaceFluxSolver
    %SWEFACEFLUXSOLVER2D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ fluxM ] = evaluate( obj, hmin, gra, nx, ny, fm)
            [ fluxM ] = mxEvaluateSurfFlux( hmin, gra, nx, ny, fm);
        end
    end
    
end

