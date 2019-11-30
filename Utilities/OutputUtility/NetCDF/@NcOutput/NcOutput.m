classdef NcOutput < AbstractNcOutput
    
    properties ( SetAccess = protected )
        %> output NetCDF file
        ncfile
        timeVarableId
        fieldVarableId
        filename
        vtkOutput
    end
    
    methods
        function obj = NcOutput( mesh, casename, Nfield, dt, varIndex )
            obj = obj@AbstractNcOutput( mesh, casename, Nfield, dt, varIndex );
        end
        
        %> create NetCDF output file
        initFromMesh( obj, mesh, fileName, outputIntervalNum, varIndex );
        %> output result
        outputResult( obj, time, field );
        
        [ field ] = readOutputResult( obj, step );
        
        writeResultToVtk( obj, step, field );
        writeOutputResultToVtk( obj, step );
        
        function closeOutputFile( obj )
            obj.ncfile.delete();
%             obj.outputTime = ncread( obj.filename, 'time' );
        end
    end
    
end

