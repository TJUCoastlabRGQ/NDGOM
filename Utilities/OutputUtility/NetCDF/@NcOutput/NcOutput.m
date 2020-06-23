classdef NcOutput < AbstractNcOutput
    
    properties ( SetAccess = protected )
        %> output NetCDF file
        ncfile
        timeVarableId
        fieldVarableId
        vtkOutput
    end
    
    properties
        %> order of the nc file to be written
        fileOrder
        %> Step contained in each nc file
        StepPerFile
        %> Number of the nc files
        Numfile
        %> file name of NetCDF file
        filename
        %> true for file is open
        isOpen
        %> ID of ncfile
        ncid
        %> Index of the ncfile been written at present
        OnIndex
    end
    
    methods
        function obj = NcOutput( mesh, casename, Nfield, dt, varIndex )
            obj = obj@AbstractNcOutput( mesh, casename, Nfield, dt, varIndex );
        end
        
        %> create NetCDF output file
        initFromMesh( obj, mesh, filename, casename, outputIntervalNum, Nfield, VarIndex, fieldName, dt );
        %> output result
        outputResult( obj, time, field );
        
        [ field ] = readOutputResult( obj, step );
        
%         writeResultToVtk( obj, step, field );
        writeOutputResultAtStepToVtk( obj, step );
        
        writeOutputResultAtTimePointToVtk(obj, timePoint);
        
        function closeOutputFile(obj)
            obj.delete();
        end
        
        function delete( obj )
            for n = 1:numel(obj.ncid)
               if obj.isOpen(n) ~= false
                obj.isOpen(n) = false;
                netcdf.close(obj.ncid(n));
               end
            end
        end% func
        
        function outputFinalResult( obj, time, field )
            obj.outputResult( time, field);
            if obj.isOpen(obj.fileOrder) == true
                obj.isOpen(obj.fileOrder) = false;
                netcdf.close(obj.ncid(obj.fileOrder));
            end
        end
        
        function closeFile( obj, varargin  )
            if(nargin == 0)
                obj.isOpen(obj.fileOrder) = false;
                netcdf.close( obj.ncid(obj.fileOrder) );
            elseif(obj.isOpen(int64(varargin{1})) == true)
                obj.isOpen(int64(varargin{1})) = false;
                netcdf.close( obj.ncid(varargin{1} ));
            end
        end
        
        
        %> The merge output result is to be added here
        function mergeOutputResult(obj)
            %> open netcdf file first
            for i = 1:numel(obj.filename)
                if obj.isOpen(i) == false
                    obj.ncid(i) = netcdf.open( obj.filename{i}, 'WRITE');
                    obj.isOpen(i) = true;
                end
            end
            obj.outputStep = numel(netcdf.getVar(obj.ncid(1),0)); % Get number of time points
            for i = 2:numel(obj.filename)
                Time = netcdf.getVar(obj.ncid(i),0);
                %                 field = netcdf.getVar(obj.ncfile.ncid(i),1);
                Info = ncinfo(obj.filename{i});
                for n = 1:numel(Time)
                    startInd = obj.outputStep;
                    countInd = 1;
                    netcdf.putVar(obj.ncid(1), obj.timeVarableId, startInd, countInd, Time(n));
                    for m = 1:numel(Info.Variables) - 1
                        field = netcdf.getVar(obj.ncid(i),m);
                        % output physical field
                        startInd = [ 0, 0, 0, obj.outputStep ];
                        countInd = [ Info.Variables( m+1 ).Size(1),  Info.Variables( m+1 ).Size(2),  Info.Variables( m+1 ).Size(3), 1 ];
                        netcdf.putVar(obj.ncid(1), obj.fieldVarableId(m), startInd, countInd, field(:,:,:,n));
                        % increase output step num
                        
                    end
                    obj.outputStep = obj.outputStep + 1;
                end
                obj.isOpen(i) = false;
                netcdf.close(obj.ncid(i));
                str = obj.filename(i);
                delete(str{1});
            end
            obj.closeFile(1);
        end
        
        function defineIntoNetcdfFile( obj, index )
            obj.ncid(index) = netcdf.create( obj.filename{index}, 'CLOBBER');
            obj.isOpen(index) = true;
            for n = 1:numel(obj.ncfile.ncDim)
                obj.ncfile.ncDim(n).defineIntoNetcdfFile( obj.ncid(index) );
            end
            
            for n = 1:numel(obj.ncfile.ncVar)
                obj.ncfile.ncVar(n).defineIntoNetcdfFile( obj.ncid(index) );
            end
            netcdf.endDef(obj.ncid(index));
        end
        
    end
    
end

