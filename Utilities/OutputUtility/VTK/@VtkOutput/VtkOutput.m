classdef VtkOutput < AbstractVtkOutput
    
    methods (Access = public)
        function obj = VtkOutput( mesh, fieldName, casename, Nfield, dt, varIndex )
            obj = obj@AbstractVtkOutput( mesh, fieldName, casename, Nfield, dt, varIndex );
        end
        
        initFromMesh( obj, mesh );
        
        
        
        function closeOutputFile(obj)
            obj.delete();
        end
        %For vtk, doing nothing
        
        
        %> output
        function outputIntervalResult( obj, time, field2d )
            if ( time - obj.timePrevious ) > obj.timeInterval
                obj.outputResult( field2d );
                obj.timePrevious = time;
            end
        end
        
        function outputFinalResult( obj, ~, field2d )
            obj.outputResult( field2d );
        end
        
        % read output
        readOutputResult( obj, timeStep )
    end
end