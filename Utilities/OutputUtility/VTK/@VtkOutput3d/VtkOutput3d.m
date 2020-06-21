classdef VtkOutput3d < VtkOutput2d
    properties
        mesh3d
        varIndex3d
        
        Nfield3d
        varName3d
        outputStep3d
        timePrevious3d
    end
    
    properties
        Np3d
        SEToV3d
        CellVertList3d
        Npoint3d
        Points3d
        Ncell3d
        Ncon3d
        ctype3d
    end
    
    methods (Access = public)
        function obj = VtkOutput3d( physMat, mesh2d, mesh3d, casename2d,  Nfield2d, Nfield3d, dt, varIndex2d, varIndex3d )
            obj = obj@VtkOutput2d( physMat, mesh2d, casename2d, Nfield2d, dt, varIndex2d );
            obj.mesh3d = mesh3d;
            obj.varIndex3d = varIndex3d;
            obj.Nfield3d = Nfield3d;
            obj.outputStep3d = 0;
            obj.timePrevious3d = 0;
            obj.varName3d = cell(numel(varIndex3d),1);
            for i = 1:numel(varIndex3d)
                 obj.varName3d{i} = physMat.fieldName3d{varIndex3d(i)};
            end
        end

        initFromMesh( obj, mesh2d, mesh3d );
        % drawResult( obj )
        
        readOutputResult( obj, timeStep )
        
        %> output
        function outputIntervalResult( obj, time, field2d, field3d )
            if ( time - obj.timePrevious3d ) > obj.timeInterval
                obj.outputResult( time, field2d, field3d );
                obj.timePrevious3d = time;
            end
        end
        
        function outputStepResult( obj, field2d, field3d )
            obj.outputResult( [], field2d, field3d );
        end
        
        function outputFinalResult( obj, time, field2d, field3d )
            obj.outputResult( time, field2d, field3d);
        end        
    end
end