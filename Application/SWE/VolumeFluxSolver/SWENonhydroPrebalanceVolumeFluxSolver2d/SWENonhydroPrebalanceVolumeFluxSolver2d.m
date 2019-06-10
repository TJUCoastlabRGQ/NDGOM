classdef SWENonhydroPrebalanceVolumeFluxSolver2d < SWEPrebalanceVolumeFlux2d & ...
        SWENonhydroVolumeFluxSolver2d
    %SWENONHYDROPREBALANCEVOLUMEFLUXSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
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

