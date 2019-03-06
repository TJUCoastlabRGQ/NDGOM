classdef NcOutput < AbstractOutputFile
    
    properties ( SetAccess = protected )
        %> output NetCDF file
        ncfile
        timeVarableId
        fieldVarableId
        filename
        vtkOutput
    end
    
    methods
        function obj = NcOutput( casename, Nfield, dt )
            obj = obj@AbstractOutputFile( casename, Nfield, dt );
        end
        
        %> create NetCDF output file
        initFromMesh( obj, mesh, filename, outputIntervalNum, varIndex );
        %> output result
        outputResult( obj, time, field );
        
        [ field ] = readOutputResult( obj, step );
        
        writeResultToVtk( obj, step, field );
        writeOutputResultToVtk( obj, step );
        
        function mergeOutputResult(obj)
            obj.ncfile.openNetcdfFile;
            obj.outputStep = numel(netcdf.getVar(obj.ncfile.ncid(1),0)); % Get number of time points
            for i = 2:numel(obj.ncfile.fileName)
                Time = netcdf.getVar(obj.ncfile.ncid(i),0);
                field = netcdf.getVar(obj.ncfile.ncid(i),1);
                for n = 1:numel(Time)
                    startInd = obj.outputStep;
                    countInd = 1;
                    netcdf.putVar(obj.ncfile.ncid(1), obj.timeVarableId, startInd, countInd, Time(n));
                    
                    % output physical field
                    startInd = [ 0, 0, 0, obj.outputStep ];
                    countInd = [ size(field,1), size(field,2), size(field,3), 1 ];
                    netcdf.putVar(obj.ncfile.ncid(1), obj.fieldVarableId, startInd, countInd, field(:,:,:,n));
                    % increase output step num
                    obj.outputStep = obj.outputStep + 1;
                end
                obj.ncfile.closeNetcdfFile(i);
%                 str = obj.ncfile.fileName{i};
%                 delete(str);
            end
            obj.ncfile.closeNetcdfFile(1);
        end
        
        function closeOutputFile( obj )
            obj.ncfile.delete();
            obj.outputTime = ncread( obj.filename, 'time' );
        end
    end
    
end

