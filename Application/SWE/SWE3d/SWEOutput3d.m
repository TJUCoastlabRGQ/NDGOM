classdef SWEOutput3d < NcOutput
    %LSWEOUTPUT3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties ( SetAccess = protected )
        mesh2d
        mesh3d
    end
    
    methods
        function obj = SWEOutput3d( casename, Nfield, dt )
            obj = obj@NcOutput( casename, Nfield, dt );
        end
        
        function initFromMesh( obj, physMat, mesh2d, mesh3d , filename, outputIntervalNum)
            % define dimension
            dimTime = NdgNcDim('Nt', 0);
            dimK2 = NdgNcDim('K2d', mesh2d.K);
            dimNp2 = NdgNcDim('Np2d', mesh2d.cell.Np);
            dimNfield2 = NdgNcDim('Nfield2d', numel( physMat.outputFieldOrder2d ) );
            
            dimK3 = NdgNcDim('K3d', mesh3d.K );
            dimNp3 = NdgNcDim('Np3d', mesh3d.cell.Np );
            dimNfield3 = NdgNcDim('Nfield3d', numel( physMat.outputFieldOrder ) );
            
            % define variable
            varTime = NdgNcVar('time', dimTime, enumNcData.NC_DOUBLE );
            varField2 = NdgNcVar('fphys2d', ...
                [ dimNp2, dimK2, dimNfield2, dimTime], ...
                enumNcData.NC_DOUBLE);
            
            varField3 = NdgNcVar('fphys3d', ...
                [ dimNp3, dimK3, dimNfield3, dimTime], ...
                enumNcData.NC_DOUBLE);
            
            if ~isdir(obj.casename)
                mkdir(obj.casename);
            end
            
            obj.ncfile = NdgNcFile( filename, ...
                [dimTime, dimK2, dimNp2, dimNfield2, dimK3, dimNp3, dimNfield3], ...
                [varTime, varField2, varField3]);
            %
            if floor(outputIntervalNum/numel(obj.ncfile.fileName))<1
                error( 'Too many output nc file!' );
            else
                obj.ncfile.StepPerFile = floor(outputIntervalNum/numel(obj.ncfile.fileName));
            end
            %
            % obj.ncfile.varIndex = varIndex;
            % % init file
            for n = 1:numel(obj.ncfile.fileName)
                obj.ncfile.defineIntoNetcdfFile(n);
            end
            % % set properties
            obj.timeVarableId = varTime.id;
            obj.fieldVarableId(1) = varField2.id;
            obj.fieldVarableId(2) = varField3.id;
            obj.mesh2d = mesh2d;
            obj.mesh3d = mesh3d;
        end
        
        function outputResult( obj, physMat, time, field2d, field3d )
            if ( time - obj.timePrevious ) > obj.timeInterval
                putVarToNcFile( obj, physMat, time, field2d, field3d);
                if obj.outputStep == ( obj.ncfile.StepPerFile - 1 ) && obj.ncfile.fileOrder ~= obj.ncfile.Numfile
                    obj.ncfile.closeNetcdfFile(obj.ncfile.fileOrder);
                    obj.ncfile.fileOrder = obj.ncfile.fileOrder + 1;
                    obj.outputStep = 0;
                else
                    % increase output step num
                    obj.outputStep = obj.outputStep + 1;
                end
            end
        end
        
        function outputFinalResult(obj, physMat, time, field2d, field3d)
            putVarToNcFile( obj, physMat, time, field2d, field3d);
        end
        
        function [fphys2d, fphys3d] = readOutputResult( obj, physMat, timeStep )
            if timeStep <= obj.outputStep
                startInd = [ 1, 1, 1, timeStep ];
                countInd = [ obj.mesh2d.cell.Np, obj.mesh2d.K, numel( physMat.outputFieldOrder2d ), 1 ];
                fphys2d = ncread( obj.ncfile.fileName{1}, 'fphys2d', startInd, countInd );
                
                startInd = [ 1, 1, 1, timeStep ];
                countInd = [ obj.mesh3d.cell.Np, obj.mesh3d.K, numel( physMat.outputFieldOrder ), 1 ];
                fphys3d = ncread(obj.ncfile.fileName{1}, 'fphys3d', startInd, countInd );
            else
                error('Output time step is less than inptu step number!')
            end
        end
    end
    
end

function putVarToNcFile(obj, physMat, time, field2d, field3d)
startInd = obj.outputStep;
countInd = 1;
netcdf.putVar(obj.ncfile.ncid(obj.ncfile.fileOrder), obj.timeVarableId, startInd, countInd, time);

% output physical field
startInd = [ 0, 0, 0, obj.outputStep ];
countInd = [ obj.mesh2d.cell.Np, obj.mesh2d.K, numel( physMat.outputFieldOrder2d ), 1 ];
netcdf.putVar(obj.ncfile.ncid(obj.ncfile.fileOrder), obj.fieldVarableId(1), startInd, countInd, field2d(:,:,physMat.outputFieldOrder2d));

startInd = [ 0, 0, 0, obj.outputStep ];
countInd = [ obj.mesh3d.cell.Np, obj.mesh3d.K, numel( physMat.outputFieldOrder ), 1 ];
netcdf.putVar(obj.ncfile.ncid(obj.ncfile.fileOrder), obj.fieldVarableId(2), startInd, countInd, field3d(:,:,physMat.outputFieldOrder));
end

