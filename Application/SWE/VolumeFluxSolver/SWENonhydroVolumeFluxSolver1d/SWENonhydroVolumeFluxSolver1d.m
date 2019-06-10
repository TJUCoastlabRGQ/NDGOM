classdef SWENonhydroVolumeFluxSolver1d  < SWEVolumeFluxSolver1d
    %SWENONHYDROVOLUMEFLUXSOLVER1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
         function [ E ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E ] = evaluate@SWEVolumeFluxSolver1d( obj, hmin, gra, mesh, fphys );
            [ E(:,:,3) ] = obj.matGetNonhydroVerticalVolumeFlux( mesh.status,...
                fphys(:,:,2), fphys(:,:,5), fphys(:,:,1));
        end        
    end
    
    methods(Access = protected)
        
        function [ huw ] = matGetNonhydroVerticalVolumeFlux( obj, status, hu, hw, h)
            [ huw ] = mxGetNonhydroVerticalVolumeFlux( status, hu, hw, h);            
        end
        
    end    
    
end

