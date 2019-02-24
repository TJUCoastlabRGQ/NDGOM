function [ outputObj ] = initOutput( physMat, mesh )
%INITOUTPUTFILE Summary of this function goes here
%   Detailed explanation goes here

if physMat.option.isKey('outputTimeInterval')
    dt = physMat.getOption( 'outputTimeInterval' );
else
    error( 'Please set the output time interval option "outputTimeInterval".' );
end

if physMat.option.isKey('outputCaseName')
    casename = physMat.getOption( 'outputCaseName' );
else
    error( 'Please set the output case name option "outputCaseName".' );
end

if physMat.option.isKey('outputFieldOrder') 
    OutputFieldNum = numel(physMat.getOption( 'outputFieldOrder' ));
else
    OutputFieldNum = physMat.Nvar;
end

if physMat.option.isKey('outputNcfileNum') 
    OutputFileNum = physMat.getOption( 'outputNcfileNum' );
else
    OutputFileNum = 1;
end

if physMat.option.isKey('outputTimeInterval')
    finalTime = physMat.getOption( 'finalTime' );
    deltaTimeStep = physMat.getOption( 'outputTimeInterval' );
    outputIntervalNum = finalTime/deltaTimeStep;
end

if physMat.option.isKey('outputFieldOrder')
    varIndex = physMat.getOption('outputFieldOrder');
else
    varIndex = physMat.varFieldIndex;
end

if physMat.option.isKey('outputType')
    if ( physMat.getOption('outputType') == enumOutputFile.NetCDF )
        [ outputObj ] = initNcOutput( physMat, casename, mesh, dt, OutputFieldNum, OutputFileNum, outputIntervalNum, varIndex );
    elseif ( physMat.getOption('outputType') == enumOutputFile.VTK )
        [ outputObj ] = initVtkOutput( physMat, casename, mesh, dt );
    elseif ( physMat.getOption('outputType') == enumOutputFile.None )
        error( ['Please set the output file type "outputType" ', ...
            'as one of the following:\nNetCDF\nVTK\n'] );
    end
else% default output type NetCDF
    [ outputObj ] = initNcOutput( physMat, casename, mesh, dt, OutputFieldNum, OutputFileNum, outputIntervalNum, varIndex );
end

end


function [ outputObj ] = initNcOutput( physMat, casename, mesh, dt, OutputFieldNum, OutputFileNum, outputIntervalNum, varIndex ) 
outputObj = [];
for m = 1:physMat.Nmesh
    filename = [];
    for n = 1:OutputFileNum
    filename{n} = [ casename,'/',casename, '.', num2str(m), '-', num2str(physMat.Nmesh),'.', num2str(n),'.','nc' ];
    end
    outputObj = [ outputObj, NcOutput( casename, OutputFieldNum, dt ) ];
    outputObj(m).initFromMesh( mesh(m), filename, outputIntervalNum, varIndex );
end
end

function [ outputObj ] = initVtkOutput( physMat, casename, mesh, dt )
outputObj = [];

if mesh(1).type == enumMeshDim.Two
    for m = 1:physMat.Nmesh
        filename = [ casename, '.', num2str(m), '-', num2str(physMat.Nmesh) ];
        outputObj = [ outputObj, VtkOutput2d( filename, physMat.Nvar, dt ) ];
        outputObj(m).initFromMesh( mesh(m) );
    end
elseif (mesh(1).type == enumMeshDim.Three)
    for m = 1:physMat.Nmesh
        filename = [ casename, '.', num2str(m), '-', num2str(physMat.Nmesh) ];
        outputObj = [ outputObj, VtkOutput3d( filename, physMat.Nvar, dt ) ];
        outputObj(m).initFromMesh( mesh(m) );
    end
end

end

