classdef SWEVolumeFluxSolver1d < SWEAbstractVolumeFluxSolver
    %SWEVOLUMEFLUXSOLVER1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ E ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E ] = mxEvaluateFlux1d( hmin, gra, mesh.status, fphys );
        end                
    end
    
end

