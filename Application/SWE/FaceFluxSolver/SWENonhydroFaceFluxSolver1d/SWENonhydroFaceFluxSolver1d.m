classdef SWENonhydroFaceFluxSolver1d < SWEFaceFluxSolver1d
    %SWENONHYDROFACEFLUXSOLVER1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function [ fluxM ] = evaluate( obj, hmin, gra, nx, fm, mesh, edge)
            [ fluxM ] = evaluate@SWEFaceFluxSolver1d( obj, hmin, gra, nx, fm, mesh, edge);
            fluxM(:,:,3) = fm(:,:,2) .* fm(:,:,5) ./ fm(:,:,1) .* nx; 
%             fluxM(:,:,3) = mxEvaluateVerticalSurfFlux( int8(mesh.status), edge.FToE, hmin,...
%                 fm(:,:,1), fm(:,:,2), fm(:,:,5), nx );          
        end          
    end
    
end

