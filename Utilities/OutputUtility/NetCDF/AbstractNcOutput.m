classdef AbstractNcOutput < AbstractOutputFile
    %ABSTRACTNCOUTPUT 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
    end
    
    methods
        function obj = AbstractNcOutput( mesh, casename, Nfield, dt, varIndex )
            obj = obj@AbstractOutputFile( mesh, casename, Nfield, dt, varIndex );
        end
    end
    
    methods(Abstract)
        
        initFromMesh( obj, mesh, fileName, outputIntervalNum, varIndex );
                
        outputResult( obj, time, field );
        
    end
    
end

