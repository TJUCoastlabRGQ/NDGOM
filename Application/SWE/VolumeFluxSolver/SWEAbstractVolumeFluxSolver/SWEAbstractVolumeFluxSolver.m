classdef SWEAbstractVolumeFluxSolver
    %SWEABSTRACTVOLUMEFLUXSOLVER 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods(Abstract)
        [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
    end
    
end

