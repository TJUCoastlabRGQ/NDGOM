classdef SWEAbstractFaceFluxSolver
    %SWEABSTRACTFACEFLUXSOLVER �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    methods(Abstract)
        [ fluxS ] = evaluate( obj, hmin, gra, nx, ny, fm)
    end
    
end

