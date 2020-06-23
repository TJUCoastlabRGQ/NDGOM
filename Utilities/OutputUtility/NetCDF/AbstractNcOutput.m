classdef AbstractNcOutput < AbstractOutputFile
    %ABSTRACTNCOUTPUT �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
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

