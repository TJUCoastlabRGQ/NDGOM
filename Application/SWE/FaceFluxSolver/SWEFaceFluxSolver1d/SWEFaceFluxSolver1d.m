classdef SWEFaceFluxSolver1d < SWEAbstractFaceFluxSolver
    %SWEFACEFLUXSOLVER1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, fm, ~ , ~)
            [ fluxS ] = mxEvaluateSurfFlux1d( hmin, gra, nx, fm);
        end
    end
    
end

