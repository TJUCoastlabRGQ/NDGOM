classdef SWEPrebalanceVolumeFlux1d < SWEAbstractVolumeFluxSolver
    %SWEPREBALANCEVOLUMEFLUX1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    methods
        function [ E ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E ] = mxEvaluateFlux1d( hmin, gra, mesh.status, fphys );
        end             
    end
    
end

