classdef SWEFaceFluxSolver1d < SWEAbstractFaceFluxSolver
    %SWEFACEFLUXSOLVER1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, fm, ~ , ~)
            [ fluxS ] = mxEvaluateSurfFlux1d( hmin, gra, nx, fm);
        end
    end
    
end

