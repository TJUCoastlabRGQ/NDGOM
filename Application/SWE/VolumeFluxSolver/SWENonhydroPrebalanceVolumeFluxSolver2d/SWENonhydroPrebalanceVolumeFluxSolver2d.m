classdef SWENonhydroPrebalanceVolumeFluxSolver2d < SWEPrebalanceVolumeFlux
    %SWENONHYDROPREBALANCEVOLUMEFLUXSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
            [ E, G ] = evaluate@SWEPrebalanceVolumeFlux( obj, hmin, gra, mesh.status, fphys );
            E(:,:,4) = fphys(:,:,2) .* fphys(:,:,6) ./ fphys(:,:,1);
            G(:,:,4) = fphys(:,:,3) .* fphys(:,:,6) ./ fphys(:,:,1);                
        end         
    end
    
end

