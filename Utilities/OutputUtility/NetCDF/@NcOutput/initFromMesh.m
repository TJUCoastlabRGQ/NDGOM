function initFromMesh( obj, filename, outputIntervalNum, varIndex )

% set vtk output
% if (obj.mesh.type == enumMeshDim.One)
%     obj.vtkOutput = VtkOutput1d(obj.casename, obj.Nfield, obj.timeInterval);
% elseif (obj.mesh.type == enumMeshDim.Two)
%     obj.vtkOutput = VtkOutput2d(obj.mesh, obj.casename, obj.Nfield, obj.timeInterval, varIndex);
% elseif (obj.mesh.type == enumMeshDim.Three)
%     obj.vtkOutput = VtkOutput3d(obj.mesh3d, obj.casename, obj.Nfield3d, obj.timeInterval3d);
% end

% obj.vtkOutput.initFromMesh( mesh );

% define dimension
dimTime = NdgNcDim('Nt', 0);
dimK = NdgNcDim('K', obj.mesh.K);
dimNp = NdgNcDim('Np', obj.mesh.cell.Np);
dimNfield = NdgNcDim('Nvar', obj.Nfield);
% dimNfield = NdgNcDim('Nvar', 6);
% define variable
varTime = NdgNcVar('time', dimTime, enumNcData.NC_DOUBLE );
varField = NdgNcVar('fphys', [dimNp, dimK, dimNfield, dimTime], enumNcData.NC_DOUBLE);


obj.ncid = zeros(size(filename));
obj.isOpen = false * ones(size(filename));
obj.fileOrder = 1;
obj.Numfile = numel(filename);
obj.fileName = filename;
% define file
% obj.filename = [ obj.casename, '/', obj.casename, '.nc' ];
obj.ncfile = NdgNcFile( obj, ...
    [dimTime, dimK, dimNp, dimNfield], [varTime, varField]);

if floor(outputIntervalNum/numel(obj.fileName))<1
    error( 'Too many output nc file!' );
else
    obj.StepPerFile = floor(outputIntervalNum/numel(obj.fileName));
end

% obj.ncfile.varIndex = varIndex;
% init file
for n = 1:numel(obj.fileName)
obj.defineIntoNetcdfFile(n);
end

% set properties
obj.timeVarableId = varTime.id;
obj.fieldVarableId = varField.id;
% obj.mesh = mesh;

end