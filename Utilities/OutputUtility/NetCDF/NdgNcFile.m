classdef NdgNcFile < handle
    
    properties(SetAccess = protected)
        %         %> true for file is open
        %         isOpen
        %         %> ID of ncfile
        %         ncid
        %> array of dimensions in NetCDF file
        ncDim
        %> array of variables in NetCDF file
        ncVar
        
        NcOutPut
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
        function obj = NdgNcFile( NcOutPut, ncdim, ncvar )
            obj.NcOutPut = NcOutPut;
            obj.ncDim = ncdim;
            obj.ncVar = ncvar;
            %             obj.ncid = zeros(size(filename));
            %             obj.isOpen = false * ones(size(filename));
            %             obj.fileOrder = 1;
            %             obj.Numfile = numel(filename);
        end% func
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
