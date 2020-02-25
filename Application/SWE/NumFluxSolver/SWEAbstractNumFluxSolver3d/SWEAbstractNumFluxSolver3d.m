%> Abstract solver for evaluating numerical flux term
classdef SWEAbstractNumFluxSolver3d
    
    methods(Abstract)
        [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, fm, fp );
    end
    
end
