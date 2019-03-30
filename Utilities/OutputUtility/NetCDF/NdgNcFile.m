classdef NdgNcFile < handle
    
    properties(SetAccess = protected)
        %> true for file is open
        isOpen
        %> ID of ncfile
        ncid
        %> array of dimensions in NetCDF file
        ncDim
        %> array of variables in NetCDF file
        ncVar
        %> file name of NetCDF file
        fileName
    end
    
    properties
        %> order of the nc file to be written
        fileOrder
        %> Step contained in each nc file
        StepPerFile
        %> Index of variable to be output
        varIndex
        %> Number of the nc files
        Numfile
    end
    
    methods
        %======================================================================
        %> @brief Brief description of the function
        %>
        %> More detailed description.
        %>
        %> @param arg1 First argument
        %> @param arg2 Second argument
        %>
        %> @retval out1 return value for the first output variable
        %> @retval out2 return value for the second output variable
        %======================================================================
        %> This function is part of the NDGOM software.
        %> @author li12242, Tianjin University, li12242@tju.edu.cn
        %======================================================================
        function obj = NdgNcFile( filename, ncdim, ncvar )
            obj.fileName = filename;
            obj.ncDim = ncdim;
            obj.ncVar = ncvar;
            obj.ncid = zeros(size(filename));
            obj.isOpen = false * ones(size(filename));
            obj.fileOrder = 1;
            obj.Numfile = numel(filename);
        end% func
        
        function delete( obj )
            for n = 1:numel(obj.ncid)
                obj.closeNetcdfFile(n);
            end
        end% func
        
        function closeNetcdfFile( obj, index )
            if(obj.isOpen(index)) % if netcdf file is still open
                obj.isOpen(index) = false;
                netcdf.close( obj.ncid(index) );
            end
        end
        
        
        function openNetcdfFile(obj)
            for i = 1:numel(obj.fileName)
                if obj.isOpen(i) == false
                    obj.ncid(i) = netcdf.open( obj.fileName{i}, 'WRITE');
                    obj.isOpen(i) = true;
                end
            end
        end
        
        function defineIntoNetcdfFile( obj, index )
            obj.ncid(index) = netcdf.create( obj.fileName{index}, 'CLOBBER');
            obj.isOpen(index) = true;
            for n = 1:numel(obj.ncDim)
                obj.ncDim(n).defineIntoNetcdfFile( obj.ncid(index) );
            end
            
            for n = 1:numel(obj.ncVar)
                obj.ncVar(n).defineIntoNetcdfFile( obj.ncid(index) );
            end
            netcdf.endDef(obj.ncid(index));
        end
        
        function deleteNetcdfFile(obj, index)
            str = obj.fileName(index);
            delete(str{1});
        end
    end
end

function [ncid, ncDim, ncVar] = readFromNetcdfFile( fileName )
ncid = netcdf.open(fileName, 'NOWRITE');

dimId = netcdf.inqDimIDs( ncid );
Ndim = numel(dimId);
ncdim = cell(Ndim, 1);
for n = 1:Ndim
    % read dimension name and ID
    [ dimName, dimLength ] = netcdf.inqDim( ncid, dimId(n) );
    dim = NdgNcDim( dimName, dimLength );
    dim.setDimId( dimId(n) );
    ncdim{n} = dim;
end
ncDim = [ncdim{:}];

varId = netcdf.inqVarIDs( ncid );
Nvar = numel(varId);
ncvar = cell(Nvar, 1);
for n = 1:Nvar
    % read variable name and ID
    [varName, dataType, dimIds, ~] = netcdf.inqVar(ncid, varId(n) );
    ncdim = ncDim( dimIds + 1 );
    var = NdgNcVar( varName, ncdim, NdgNcType( dataType ) );
    var.setVarId( varId(n) );
    ncvar{n} = var;
end
ncVar = [ ncvar{:} ];

end% func
