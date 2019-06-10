classdef SWEVolumeFluxSolver1d < SWEAbstractVolumeFluxSolver
    %SWEVOLUMEFLUXSOLVER1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function [ E ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E ] = mxEvaluateFlux1d( hmin, gra, mesh.status, fphys );
        end                
    end
    
end

