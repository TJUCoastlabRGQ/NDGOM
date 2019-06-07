classdef SWENonhydroFaceFluxSolver1d < SWEFaceFluxSolver1d
    %SWENONHYDROFACEFLUXSOLVER1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ fluxM ] = evaluate( obj, hmin, gra, nx, fm)
            [ fluxM ] = evaluate@SWEFaceFluxSolver1d( obj, hmin, gra, nx, fm);
            fluxM(:,:,3) = fm(:,:,2) .* fm(:,:,5) ./ fm(:,:,1) .* nx;            
        end          
    end
    
end

