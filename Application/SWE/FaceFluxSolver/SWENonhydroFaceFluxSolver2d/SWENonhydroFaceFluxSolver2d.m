classdef SWENonhydroFaceFluxSolver2d < SWEFaceFluxSolver2d
    %SWENONHYDROFACEFLUXSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function [ fluxM ] = evaluate( obj, hmin, gra, nx, ny, fm)
            [ fluxM ] = evaluate@SWEFaceFluxSolver2d( obj, hmin, gra, nx, ny, fm);
            fluxM(:,:,4) = fm(:,:,2) .* fm(:,:,6) ./ fm(:,:,1) .* nx + ...
                fm(:,:,3) .* fm(:,:,6) ./ fm(:,:,1) .* ny;            
        end        
    end
    
end

