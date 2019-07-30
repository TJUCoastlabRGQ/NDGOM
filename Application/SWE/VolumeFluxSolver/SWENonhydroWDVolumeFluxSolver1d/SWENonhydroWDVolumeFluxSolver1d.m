classdef SWENonhydroWDVolumeFluxSolver1d < SWEWDVolumeFluxSolver1d & ...
        SWENonhydroVolumeFluxSolver1d
    %SWENONHYDROPREBALANCEVOLUMEFLUXSOLVER1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ E ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E ] = evaluate@SWEWDVolumeFluxSolver1d( obj, hmin, gra, mesh, fphys );
            [ E(:,:,3) ] = obj.matGetNonhydroVerticalVolumeFlux(...
                mesh.status, fphys(:,:,2), fphys(:,:,5), fphys(:,:,1));        
        end             
    end
    
end
