classdef Visual1d < AbstractVisual
    %VISUAL1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
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

