function [ outputObj ] = matInitOutput( obj, mesh, fieldName )
%INITOUTPUTFILE Summary of this function goes here
%   Detailed explanation goes here
%> At present, multiple-file output is only supported for outputIntervalType with value
%> enumOutputInterval.DeltaTime.  For enumOutputInterval.DeltaStep, this is not surpported
%> since the number of files and the output time steps for each file are not known in advance
if obj.getOption('outputIntervalType') == enumOutputInterval.DeltaTime
    if obj.option.isKey('outputTimeInterval')
        dt = obj.getOption( 'outputTimeInterval' );
        tstart = obj.getOption('startTime');
        tend = obj.getOption('finalTime');
        outputIntervalNum = ceil((tend - tstart)/dt);
    else
        error([ 'The outputIntervalType is set to be enumOutputInterval.DeltaTime,'...
            'please set the output time interval option "outputTimeInterval".'] );
    end
end

if obj.option.isKey('outputCaseName')
    casename = obj.getOption( 'outputCaseName' );
else
    error( 'Please set the output case name option "outputCaseName".' );
end

if ~isdir(casename)
    mkdir(casename);
end

%> Number of the output file, this is useful when we need to conduct the hot
%> start for computation demanding cases, such as non-hydrostatic case
if obj.option.isKey('outputNcfileNum')
    OutputFileNum = obj.getOption( 'outputNcfileNum' );
else
    OutputFileNum = 1;
end

if mesh.type == enumMeshDim.Three
    %     OutputFieldNum2d = numel( physMat.outputFieldOrder2d);
    %     varIndex2d = physMat.outputFieldOrder2d;
    OutputFieldNum = numel( obj.outputFieldOrder3d );
    varIndex = obj.outputFieldOrder3d;
    if obj.option.isKey('outputType')
        if ( obj.getOption('outputType') == enumOutputFile.NetCDF )
            [ outputObj ] = initNcOutput( obj, mesh, casename, OutputFieldNum, ...
                dt, OutputFileNum, outputIntervalNum, varIndex, fieldName );
        elseif ( obj.getOption('outputType') == enumOutputFile.VTK )
            [ outputObj ] = initVtkOutput( obj, mesh, casename, OutputFieldNum, dt, ...
                varIndex, fieldName );
        elseif ( obj.getOption('outputType') == enumOutputFile.None )
            error( ['Please set the output file type "outputType" ', ...
                'as one of the following:\nNetCDF\nVTK\n'] );
        end
    else% default output type NetCDF
        [ outputObj ] = initNcOutput( obj, mesh, casename, OutputFieldNum, ...
            dt, OutputFileNum, outputIntervalNum, varIndex, fieldName );
    end
elseif mesh.type == enumMeshDim.Two
    OutputFieldNum = numel( obj.outputFieldOrder2d );
    varIndex = obj.outputFieldOrder2d;
    if obj.option.isKey('outputType')
        if ( obj.getOption('outputType') == enumOutputFile.NetCDF )
            [ outputObj ] = initNcOutput( obj, mesh, casename, OutputFieldNum, ...
                dt, OutputFileNum, outputIntervalNum, varIndex, fieldName );
        elseif ( obj.getOption('outputType') == enumOutputFile.VTK )
            [ outputObj ] = initVtkOutput( obj, mesh, casename, OutputFieldNum, dt, ...
                varIndex, fieldName );
        elseif ( obj.getOption('outputType') == enumOutputFile.None )
            error( ['Please set the output file type "outputType" ', ...
                'as one of the following:\nNetCDF\nVTK\n'] );
        end
    else% default output type NetCDF
        [ outputObj ] = initNcOutput( obj, mesh, casename, OutputFieldNum, ...
            dt, OutputFileNum, outputIntervalNum, varIndex, fieldName );
    end
else
    %One dimension problem is not considered at present
end

end


function [ outputObj ] = initNcOutput( obj, mesh, casename, OutputFieldNum, dt, OutputFileNum, outputIntervalNum, varIndex, fieldName )
outputObj = [];
for m = 1:obj.Nmesh
    filename = cell(OutputFileNum, 1);
    if mesh.type == enumMeshDim.Three
        str = '/3d/';
    elseif mesh.type == enumMeshDim.Two
        str = '/2d/';
    end
    for n = 1:OutputFileNum
        filename{n} = [ casename, str ,casename, '.', num2str(m), '-', num2str(obj.Nmesh),'.', num2str(n),'.','nc' ];
    end
    outputObj = [ outputObj, NcOutput( mesh, casename, OutputFieldNum, dt, varIndex ) ];
    outputObj(m).initFromMesh( mesh, filename, casename, outputIntervalNum, OutputFieldNum, varIndex, fieldName, dt );
end
end

function [ outputObj ] = initVtkOutput( obj, mesh, casename, OutputFieldNum, dt, varIndex, fieldName )
outputObj = [];
if mesh.type == enumMeshDim.Three
    if ~isdir([casename,'/3d'])
        mkdir([casename,'/3d']);
    end
elseif mesh.type == enumMeshDim.Two
    if ~isdir([casename,'/2d'])
        mkdir([casename,'/2d']);
    end
end
for m = 1:obj.Nmesh
    outputObj = [ outputObj, VtkOutput( mesh, fieldName, casename,...
        OutputFieldNum,  dt, varIndex ) ];
    outputObj(m).initFromMesh( mesh );
end
end

