classdef SWEHLLCNumFluxSolver2d < SWEAbstractNumFluxSolver2d
    
    methods
        function [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, fm, fp, Nvar, ~, ~ )
            fluxS = mxEvaluate( hmin, gra, nx, ny, fm, fp );
%             fluxS = mxEvaluateNew( hmin, gra, nx, ny, fm, fp, Nvar );
%             fluxS = mxEvaluateNew( hmin, gra, nx, ny, fm, fp, 3 );
        end
    end
    
end

