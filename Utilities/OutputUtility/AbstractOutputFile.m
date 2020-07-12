classdef AbstractOutputFile < handle
    %OUTPUTABSTRACT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( SetAccess = protected )
        mesh
        Nfield
        casename
        %         outputTime
        varIndex     
    end
    
    properties
        timeInterval
        timePrevious
        outputStep
    end
    
    methods
        function obj = AbstractOutputFile( mesh, casename, Nfield, timeInterval, varIndex )
            obj.mesh = mesh;
            obj.Nfield = Nfield;
            obj.casename = casename;
            obj.timeInterval = timeInterval;
            obj.timePrevious = 0;
            obj.outputStep = 0;
            obj.varIndex = varIndex;
        end
        
        %> output
        function outputIntervalResult( obj, time, field )
            if ( time - obj.timePrevious ) > obj.timeInterval
                obj.outputResult( time, field );
                obj.timePrevious = time;
            end
        end
        
        function outputFinalResult( obj, time, field )
            obj.outputResult( time, field );
        end
    end
    
    methods ( Abstract )
        initFromMesh( obj, mesh, timeInterval )
        outputResult( obj, time, field )
        closeOutputFile( obj )
        readOutputResult( obj, timeStep )
    end
    
end

