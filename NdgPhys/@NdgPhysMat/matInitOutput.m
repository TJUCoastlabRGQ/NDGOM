function [ outputObj ] = matInitOutput( physMat )
%INITOUTPUTFILE Summary of this function goes here
%   Detailed explanation goes here
%> At present, multiple-file output is only supported for outputIntervalType with value
%> enumOutputInterval.DeltaTime.  For enumOutputInterval.DeltaStep, this is not surpported
%> since the number of files and the output time steps for each file are not known in advance
if physMat.getOption('outputIntervalType') == enumOutputInterval.DeltaTime
    if physMat.option.isKey('outputTimeInterval')
        dt = physMat.getOption( 'outputTimeInterval' );
        tstart = physMat.getOption('startTime');
        tend = physMat.getOption('finalTime');
        outputIntervalNum = ceil((tend - tstart)/dt);
    else
        error([ 'The outputIntervalType is set to be enumOutputInterval.DeltaTime,'...
            'please set the output time interval option "outputTimeInterval".'] );
    end
end

if physMat.option.isKey('outputCaseName')
    casename = physMat.getOption( 'outputCaseName' );
else
    error( 'Please set the output case name option "outputCaseName".' );
end

if ~isdir(casename)
    mkdir(casename);
end

%> Number of the output file, this is useful when we need to conduct the hot
%> start for computation demanding cases, such as non-hydrostatic case
if physMat.option.isKey('outputNcfileNum')
    OutputFileNum = physMat.getOption( 'outputNcfileNum' );
else
    OutputFileNum = 1;
end

if physMat.meshUnion(1).type == enumMeshDim.Three
    OutputFieldNum2d = numel( physMat.outputFieldOrder2d);
    varIndex2d = physMat.outputFieldOrder2d;
    OutputFieldNum3d = numel( physMat.outputFieldOrder );
    varIndex3d = physMat.outputFieldOrder;
    if physMat.option.isKey('outputType')
        if ( physMat.getOption('outputType') == enumOutputFile.NetCDF )
            [ outputObj ] = initNcOutput3d( physMat, casename, dt, OutputFieldNum2d, ...
                OutputFieldNum3d, OutputFileNum, outputIntervalNum, varIndex2d, varIndex3d );
        elseif ( physMat.getOption('outputType') == enumOutputFile.VTK )
            [ outputObj ] = initVtkOutput3d( physMat, casename, dt, OutputFieldNum2d, ...
                OutputFieldNum3d, varIndex2d, varIndex3d );
        elseif ( physMat.getOption('outputType') == enumOutputFile.None )
            error( ['Please set the output file type "outputType" ', ...
                'as one of the following:\nNetCDF\nVTK\n'] );
        end
    else% default output type NetCDF
        [ outputObj ] = initNcOutput3d( physMat, casename, dt, OutputFieldNum2d, ...
            OutputFieldNum3d, OutputFileNum, outputIntervalNum, varIndex2d, varIndex3d );
    end
else
    OutputFieldNum = numel( physMat.outputFieldOrder );
    varIndex = physMat.outputFieldOrder;
    if physMat.option.isKey('outputType')
        if ( physMat.getOption('outputType') == enumOutputFile.NetCDF )
            [ outputObj ] = initNcOutput( physMat, casename, dt, OutputFieldNum, OutputFileNum, outputIntervalNum, varIndex );
            %> we point out that, VTK file format is only valid for two-dimensional and
            %> three-dimensional cases
        elseif ( physMat.getOption('outputType') == enumOutputFile.VTK )
            [ outputObj ] = initVtkOutput2d( physMat, casename, dt, OutputFieldNum, varIndex );
        elseif ( physMat.getOption('outputType') == enumOutputFile.None )
            error( ['Please set the output file type "outputType" ', ...
                'as one of the following:\nNetCDF\nVTK\n'] );
        end
    else% default output type NetCDF
        [ outputObj ] = initNcOutput( physMat, casename, dt, OutputFieldNum, OutputFileNum, outputIntervalNum, varIndex );
    end
end

end


function [ outputObj ] = initNcOutput( physMat, casename, dt, OutputFieldNum, OutputFileNum, outputIntervalNum, varIndex )
outputObj = [];
for m = 1:physMat.Nmesh
    filename = cell(OutputFileNum, 1);
    for n = 1:OutputFileNum
        filename{n} = [ casename,'/',casename, '.', num2str(m), '-', num2str(physMat.Nmesh),'.', num2str(n),'.','nc' ];
    end
    outputObj = [ outputObj, NcOutput( physMat.meshUnion(m), casename, OutputFieldNum, dt, varIndex ) ];
    outputObj(m).initFromMesh( filename, outputIntervalNum, varIndex );
end
end

function [ outputObj ] = initNcOutput3d( physMat, casename, dt, OutputFieldNum2d, OutputFieldNum3d, OutputFileNum, outputIntervalNum, varIndex2d, varIndex3d )
outputObj = [];
if ~isdir([casename,'/2d'])
    mkdir([casename,'/2d']);
end
if ~isdir([casename,'/3d'])
    mkdir([casename,'/3d']);
end
for m = 1:physMat.Nmesh
    filename2d = cell(OutputFileNum, 1);
    filename3d = cell(OutputFileNum, 1);
    for n = 1:OutputFileNum
        filename2d{n} = [ casename,'/2d/',casename,'-2d', '.', num2str(m), '-', num2str(physMat.Nmesh),'.', num2str(n),'.','nc' ];
        filename3d{n} = [ casename,'/3d/',casename,'-3d', '.', num2str(m), '-', num2str(physMat.Nmesh),'.', num2str(n),'.','nc' ];
    end
    outputObj = [ outputObj, NcOutput3d( physMat, casename, OutputFieldNum2d, OutputFieldNum3d, dt, varIndex2d, varIndex3d ) ];
    outputObj(m).initFromMesh( physMat, filename2d, filename3d, outputIntervalNum, varIndex2d, varIndex3d );
end
end

function [ outputObj ] = initVtkOutput2d( physMat, casename, dt, OutputFieldNum, varIndex )
outputObj = [];
if ~isdir([casename,'/2d'])
    mkdir([casename,'/2d']);
end
% casename2d = [casename,'-2d'];
for m = 1:physMat.Nmesh
    outputObj = [ outputObj, VtkOutput2d( physMat, physMat.meshUnion(m), casename, OutputFieldNum, dt, varIndex ) ];
    outputObj(m).initFromMesh( physMat.meshUnion(m) );
end

end

function [ outputObj ] = initVtkOutput3d( physMat, casename, dt, OutputFieldNum2d, OutputFieldNum3d, varIndex2d, varIndex3d )
outputObj = [];
if ~isdir([casename,'/2d'])
    mkdir([casename,'/2d']);
end
if ~isdir([casename,'/3d'])
    mkdir([casename,'/3d']);
end
% casename2d = [casename,'-2d'];
% casename3d = [casename,'-3d'];
for m = 1:physMat.Nmesh
    outputObj = [ outputObj, VtkOutput3d( physMat, physMat.mesh2d(m), physMat.meshUnion(m), casename,...
        OutputFieldNum2d, OutputFieldNum3d, dt, varIndex2d, varIndex3d ) ];
    outputObj(m).initFromMesh( physMat.mesh2d(m), physMat.meshUnion(m) );
end
end

