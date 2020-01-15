classdef NcOutput < AbstractNcOutput
    
    properties ( SetAccess = protected )
        %> output NetCDF file
        ncfile
        timeVarableId
        fieldVarableId
        filename
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
        fileName
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
        initFromMesh( obj, physMat, mesh, fileName, outputIntervalNum, varIndex );
        %> output result
        outputResult( obj, time, field );
        
        [ field ] = readOutputResult( obj, step );
        
        writeResultToVtk( obj, step, field );
        writeOutputResultToVtk( obj, step );
        
        function closeOutputFile( obj )
            obj.delete();
        end
        
        function delete( obj )
            for n = 1:numel(obj.ncid)
                obj.isOpen(n) = false;
                netcdf.close(obj.ncid(n));
            end
        end% func
        
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
            for i = 1:numel(obj.fileName)
                if obj.isOpen(i) == false
                    obj.ncid(i) = netcdf.open( obj.fileName{i}, 'WRITE');
                    obj.isOpen(i) = true;
                end
            end
            obj.outputStep = numel(netcdf.getVar(obj.ncid(1),0)); % Get number of time points
            for i = 2:numel(obj.fileName)
                Time = netcdf.getVar(obj.ncid(i),0);
                %                 field = netcdf.getVar(obj.ncfile.ncid(i),1);
                Info = ncinfo(obj.fileName{i});
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
                str = obj.fileName(i);
                delete(str{1});
            end
            obj.closeFile(1);
        end
        
        function defineIntoNetcdfFile( obj, index )
            obj.ncid(index) = netcdf.create( obj.fileName{index}, 'CLOBBER');
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

