classdef SWENonhydroFaceFluxSolver1d < SWEFaceFluxSolver1d
    %SWENONHYDROFACEFLUXSOLVER1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
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

