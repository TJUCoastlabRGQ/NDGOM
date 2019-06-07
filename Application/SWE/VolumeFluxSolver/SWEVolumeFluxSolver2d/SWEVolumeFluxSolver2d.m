classdef SWEVolumeFluxSolver2d < SWEAbstractVolumeFluxSolver
    %SWEVOLUMEFLUXSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E, G ] = mxEvaluateFlux2d( hmin, gra, mesh.status, fphys );
        end        
    end
    
end

