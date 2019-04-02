function initOutput( obj, mesh2d, mesh3d )
%INITOUTPUT Summary of this function goes here
%   Detailed explanation goes here

% obj.outputFile = SWEOutput3d( obj.casename, obj.Nfield2d, obj.outputTimeInterval );
% obj.outputFile.initFromMesh( mesh(m), filename, outputIntervalNum, varIndex );
obj.option = obj.setOption( obj.option );

casename = obj.getOption('outputCaseName');
OutputFileNum = obj.getOption('outputNcfileNum');
outputIntervalNum = obj.getOption('finalTime')./obj.getOption('outputTimeInterval');
dt = obj.getOption('outputTimeInterval');
obj.outputFile = [];
for m = 1:obj.Nmesh
    filename = [];
    for n = 1:OutputFileNum
    filename{n} = [ casename,'/',casename, '.', num2str(m), '-', num2str(obj.Nmesh),'.', num2str(n),'.','nc' ];
    end
    obj.outputFile = [ obj.outputFile, SWEOutput3d( casename, 2, dt ) ];
    obj.outputFile(m).initFromMesh( obj, mesh2d, mesh3d , filename, outputIntervalNum);
end
end
