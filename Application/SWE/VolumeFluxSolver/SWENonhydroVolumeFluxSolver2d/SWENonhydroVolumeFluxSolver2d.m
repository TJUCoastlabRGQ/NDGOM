classdef SWENonhydroVolumeFluxSolver2d < SWEVolumeFluxSolver2d
    %SWENONHYDROVOLUMEFLUXSOLVER2D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E, G ] = evaluate@SWEVolumeFluxSolver2d( obj, hmin, gra, mesh, fphys );
            [ E(:,:,4), G(:,:,4) ] = obj.matGetNonhydroVerticalVolumeFlux( mesh.status,...
                fphys(:,:,2), fphys(:,:,3), fphys(:,:,6), fphys(:,:,1));
        end 
    end
    
    methods(Access = protected)
        
        function [ huw, hvw ] = matGetNonhydroVerticalVolumeFlux( obj, status, hu, hv, hw, h)
            [ huw, hvw ] = mxGetNonhydroVerticalVolumeFlux( status, hu, hv, hw, h);            
        end
        
    end
    
end

