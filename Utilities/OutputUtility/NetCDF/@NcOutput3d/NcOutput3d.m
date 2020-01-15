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
        vtkOutput3d
    end
    
    properties
        %> file name of NetCDF file
        fileName3d
        %> true for file is open
        isOpen3d
        %> ID of ncfile
        ncid3d
    end
    
    properties
        outputStep3d
        fileOrder3d
        closed
    end
    
    methods
        function obj = NcOutput3d( physMat, casename, OutputFieldNum2d, OutputFieldNum3d, dt, varIndex2d, varIndex3d )
            obj = obj@NcOutput(physMat.mesh2d(1), casename, OutputFieldNum2d, dt, varIndex2d);
            obj.Nfield3d = OutputFieldNum3d;
            obj.mesh3d = physMat.meshUnion(1);
            obj.varIndex3d = varIndex3d;
            obj.outputStep3d = 0;
            obj.closed = false;
        end
        
        %> create NetCDF output file
        initFromMesh( obj, physMat, filename2d, filename3d, outputIntervalNum, varIndex2d, varIndex3d )
        
        outputResult( obj, time, field2d, field3d );
        
        function closeOutputFile( obj )
            obj.delete();
        end
        
        function delete( obj )
            for n = 1:numel(obj.ncid3d)
                obj.isOpen(n) = false;
                obj.isOpen3d(n) = false;
                netcdf.close(obj.ncid(n));
                netcdf.close(obj.ncid3d(n));
            end
        end% func
        
        %> output
        function outputIntervalResult( obj, time, field2d, field3d )
            if ( time - obj.timePrevious ) > obj.timeInterval
                obj.outputResult( time, field2d, field3d );
                obj.timePrevious = time;
            end
        end
        
        function outputFinalResult( obj, time, field2d, field3d )
            obj.outputResult( time, field2d, field3d);
            if obj.isOpen3d(obj.fileOrder3d) == true
                obj.isOpen(obj.fileOrder) = false;
                obj.isOpen3d(obj.fileOrder3d) = false;
                netcdf.close(obj.ncid3d(obj.fileOrder3d));
                netcdf.close( obj.ncid(obj.fileOrder)   );
            end
        end
        
        
        %> The merge output result is to be added here
        function mergeOutputResult(obj)
            mergeOutputResult@NcOutput(obj);
            %> open netcdf file first
            for i = 1:numel(obj.fileName3d)
                if obj.isOpen3d(i) == false
                    obj.ncid3d(i) = netcdf.open( obj.fileName3d{i}, 'WRITE');
                    obj.isOpen3d(i) = true;
                end
            end
            obj.outputStep3d = numel(netcdf.getVar(obj.ncid3d(1),0)); % Get number of time points
            for i = 2:numel(obj.fileName3d)
                Time = netcdf.getVar(obj.ncid3d(i),0);
                %                 field = netcdf.getVar(obj.ncfile.ncid(i),1);
                Info = ncinfo(obj.fileName3d{i});
                for n = 1:numel(Time)
                    startInd = obj.outputStep;
                    countInd = 1;
                    netcdf.putVar(obj.ncid3d(1), obj.timeVarableId, startInd, countInd, Time(n));
                    for m = 1:numel(Info.Variables) - 1
                        field = netcdf.getVar(obj.ncid3d(i),m);
                        % output physical field
                        startInd = [ 0, 0, 0, obj.outputStep3d ];
                        countInd = [ Info.Variables( m+1 ).Size(1),  Info.Variables( m+1 ).Size(2),  Info.Variables( m+1 ).Size(3), 1 ];
                        netcdf.putVar(obj.ncid3d(1), obj.fieldVarableId3d(m), startInd, countInd, field(:,:,:,n));
                        % increase output step num
                        
                    end
                    obj.outputStep3d = obj.outputStep3d + 1;
                end
                obj.isOpen3d(i) = false;
                netcdf.close(obj.ncid3d(i));
                str = obj.fileName3d(i);
                delete(str{1});
            end
            obj.closeFile(1);
        end
        
        function defineIntoNetcdfFile( obj, index )
            defineIntoNetcdfFile@NcOutput( obj, index );
            obj.ncid3d(index) = netcdf.create( obj.fileName3d{index}, 'CLOBBER');
            obj.isOpen3d(index) = true;
            for n = 1:numel(obj.ncfile3d.ncDim)
                obj.ncfile3d.ncDim(n).defineIntoNetcdfFile( obj.ncid3d(index) );
            end
            
            for n = 1:numel(obj.ncfile3d.ncVar)
                obj.ncfile3d.ncVar(n).defineIntoNetcdfFile( obj.ncid3d(index) );
            end
            netcdf.endDef(obj.ncid3d(index));
        end
    end
    
end

