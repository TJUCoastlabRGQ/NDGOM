classdef SWEPrebalanceVolumeFlux2d < SWEAbstractVolumeFluxSolver
    %SWEPREBALANCEVOLUMEFLUX �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    methods
        function [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E, G ] = mxEvaluateFlux2d( hmin, gra, mesh.status, fphys );
        end        
    end
    
end

