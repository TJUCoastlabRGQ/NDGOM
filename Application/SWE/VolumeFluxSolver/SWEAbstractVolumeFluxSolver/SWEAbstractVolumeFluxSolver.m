classdef SWEAbstractVolumeFluxSolver
    %SWEABSTRACTVOLUMEFLUXSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods(Abstract)
        [ E, G ] = evaluate( obj, hmin, gra, mesh, fphys )
    end
    
end

