classdef SWENonhydroPrebalanceVolumeFluxSolver2d < SWEPrebalanceVolumeFlux2d & ...
        SWENonhydroVolumeFluxSolver2d
    %SWENONHYDROPREBALANCEVOLUMEFLUXSOLVER2D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E, G ] = evaluate@SWEPrebalanceVolumeFlux2d( obj, hmin, gra, mesh, fphys );
            [ E(:,:,4), G(:,:,4)] = obj.matGetNonhydroVerticalVolumeFlux(...
                mesh.status, fphys(:,:,2), fphys(:,:,3), fphys(:,:,6), fphys(:,:,1));        
        end         
    end
    
end

