classdef SWENonhydroHLLNumFluxSolver2d < SWEHLLNumFluxSolver2d
    %SWENONHYDROHLLNUMFLUXSOLVER2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, fm, fp, mesh, edge )
            fluxS = evaluate@SWEHLLNumFluxSolver2d( obj, hmin, gra, nx, ny, fm, fp, mesh );
            % upwind flux for the HW term
            fluxS(:,:,4) = mxEvaluateUpwindNumFlux( int8(mesh.status), edge.FToE, fm(:,:,1), fm(:,:,2),...
                fm(:,:,3), fm(:,:,6), fp(:,:,1), fp(:,:,2), fp(:,:,3), fp(:,:,6), nx, ny );
        end
    end
    
end

