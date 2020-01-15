classdef SWENonhydroWDPrebalanceVolumeFluxSolver2d < SWEWDPrebalanceVolumeFlux2d & ...
        SWENonhydroVolumeFluxSolver2d
    %SWENONHYDROPREBALANCEVOLUMEFLUXSOLVER2D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
%             [ E, G ] = mxEvaluateFlux2d( hmin, gra, mesh.status, fphys );
            [ E, G ] = evaluate@SWEWDPrebalanceVolumeFlux2d( obj, hmin, gra, mesh, fphys );
            [ E(:,:,4), G(:,:,4)] = obj.matGetNonhydroVerticalVolumeFlux(...
                mesh.status, fphys(:,:,2), fphys(:,:,3), fphys(:,:,6), fphys(:,:,1));        
        end         
    end
    
end

