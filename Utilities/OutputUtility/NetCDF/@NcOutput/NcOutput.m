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
        initFromMesh( obj, mesh, fileName, outputIntervalNum, varIndex );
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
%                 field = netcdf.getVar(obj.ncfile.ncid(i),1);
                Info = ncinfo(obj.ncfile.fileName{i});
                for n = 1:numel(Time)
                    startInd = obj.outputStep;
                    countInd = 1;
                    netcdf.putVar(obj.ncfile.ncid(1), obj.timeVarableId, startInd, countInd, Time(n));
                  for m = 1:numel(Info.Variables) - 1 
                    field = netcdf.getVar(obj.ncfile.ncid(i),m);
                    % output physical field
                    startInd = [ 0, 0, 0, obj.outputStep ];
                    countInd = [ Info.Variables( m+1 ).Size(1),  Info.Variables( m+1 ).Size(2),  Info.Variables( m+1 ).Size(3), 1 ];
                    netcdf.putVar(obj.ncfile.ncid(1), obj.fieldVarableId(m), startInd, countInd, field(:,:,:,n));
                    % increase output step num
                    
                  end
                  obj.outputStep = obj.outputStep + 1;
                end
                obj.ncfile.closeNetcdfFile(i);
                obj.ncfile.deleteNetcdfFile(i);
            end
            obj.ncfile.closeNetcdfFile(1);
        end
        
        function closeOutputFile( obj )
            obj.ncfile.delete();
%             obj.outputTime = ncread( obj.filename, 'time' );
        end
    end
    
end

