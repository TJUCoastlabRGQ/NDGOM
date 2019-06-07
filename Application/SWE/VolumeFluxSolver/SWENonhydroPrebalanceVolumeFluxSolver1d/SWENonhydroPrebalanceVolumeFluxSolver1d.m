classdef SWENonhydroPrebalanceVolumeFluxSolver1d < SWEPrebalanceVolumeFlux1d & ...
        SWENonhydroVolumeFluxSolver1d
    %SWENONHYDROPREBALANCEVOLUMEFLUXSOLVER1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function [ E ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E ] = evaluate@SWEPrebalanceVolumeFlux1d( obj, hmin, gra, mesh, fphys );
            [ E(:,:,3) ] = obj.matGetNonhydroVerticalVolumeFlux(...
                mesh.status, fphys(:,:,2), fphys(:,:,5), fphys(:,:,1));        
        end             
    end
    
end

