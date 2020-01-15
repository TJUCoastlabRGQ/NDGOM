classdef VtkOutput2d < AbstractVtkOutput
    
    properties
        varName2d
    end
    
methods (Access = public)
    function obj = VtkOutput2d( physMat, mesh, casename, Nfield, dt, varIndex )
        obj = obj@AbstractVtkOutput( mesh, casename, Nfield, dt, varIndex );
        obj.varName2d = cell(numel(varIndex),1);
        for i = 1:numel(varIndex)
            obj.varName2d{i} = physMat.fieldName2d{varIndex(i)};
        end
    end

    initFromMesh( obj, mesh );
    
    function closeOutputFile( obj )
    end
    
            %> output
        function outputIntervalResult( obj, time, field2d )
            if ( time - obj.timePrevious ) > obj.timeInterval
                obj.outputResult( time, field2d );
                obj.timePrevious = time;
            end
        end
        
        function outputFinalResult( obj, time, field2d )
            obj.outputResult( time, field2d);
        end
    
    % read output
    readOutputResult( obj, timeStep )
end
end