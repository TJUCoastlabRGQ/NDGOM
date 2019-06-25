classdef SWENonhydroFaceFluxSolver2d < SWEFaceFluxSolver2d
    %SWENONHYDROFACEFLUXSOLVER2D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function [ fluxM ] = evaluate( obj, hmin, gra, nx, ny, fm, mesh, edge)
            [ fluxM ] = evaluate@SWEFaceFluxSolver2d( obj, hmin, gra, nx, ny, fm, mesh, edge);
            fluxM(:,:,4) = mxEvaluateVerticalSurfFlux( int8(mesh.status), edge.FToE, hmin, fm(:,:,1), fm(:,:,2), fm(:,:,3),...
                fm(:,:,6), nx, ny);          
        end        
    end
    
end

