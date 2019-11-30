classdef NcOutput3d  < NcOutput
    %NCOUTPUT3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Nfield3d
        mesh3d
        varIndex3d
        
        ncfile3d
        timeVarableId3d
        fieldVarableId3d
        filename3d
        vtkOutput3d        
    end
    
    properties
        timeInterval3d
        timePrevious3d
        outputStep3d
    end    
    
    methods
        function obj = NcOutput3d( physMat, casename, OutputFieldNum2d, OutputFieldNum3d, dt, varIndex2d, varIndex3d )
            obj = obj@NcOutput(physMat.mesh2d(1), casename, OutputFieldNum2d, dt, varIndex2d);
            obj.Nfield3d = OutputFieldNum3d;
            obj.mesh3d = physMat.meshUnion(1);
            obj.timeInterval3d = dt;
            obj.timePrevious3d = 0;
            obj.outputStep3d = 0;
            obj.varIndex3d = varIndex3d;
        end
        
        %> create NetCDF output file
        initFromMesh( obj, filename2d, filename3d, outputIntervalNum, varIndex2d, varIndex3d )     
        
        outputResult( obj, time, field2d, field3d );
    end
    
end

