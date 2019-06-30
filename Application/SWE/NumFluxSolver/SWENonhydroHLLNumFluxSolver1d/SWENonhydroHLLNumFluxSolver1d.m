classdef SWENonhydroHLLNumFluxSolver1d < SWEHLLNumFluxSolver1d
    %SWENONHYDROHLLNUMFLUXSOLVER1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, fm, fp, mesh, edge )
            fluxS = evaluate@SWEHLLNumFluxSolver1d( obj, hmin, gra, nx, fm, fp, mesh, edge );
            % upwind flux for the HW term
            fluxS(:,:,3) = mxEvaluateUpwindNumFlux( int8(mesh.status), edge.FToE, fm(:,:,1), fm(:,:,2),...
                 fm(:,:,5), fp(:,:,1), fp(:,:,2), fp(:,:,5), nx );            
        end        
        
    end
    
end

