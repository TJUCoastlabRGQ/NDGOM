classdef Visual1d < AbstractVisual
    %VISUAL1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = Visual1d( mesh )
            obj = obj@AbstractVisual( mesh );
        end
        
        drawResult( obj, field );
        drawOutputResult( obj, outputObj, step, fieldId );
        drawMesh( obj );        
    end
    
end

