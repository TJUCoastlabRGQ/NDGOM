classdef SWEVolumeFluxSolver2d < SWEAbstractVolumeFluxSolver
    %SWEVOLUMEFLUXSOLVER2D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    methods
        function [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E, G ] = mxEvaluateFlux2d( hmin, gra, mesh.status, fphys );
        end        
    end
    
end

