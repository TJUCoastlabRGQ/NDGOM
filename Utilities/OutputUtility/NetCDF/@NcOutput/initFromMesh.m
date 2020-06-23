function initFromMesh( obj, mesh, filename, casename, outputIntervalNum, OutputFieldNum, VarIndex, fieldname, dt )
% ( mesh, filename, casename, outputIntervalNum, OutputFieldNum, varIndex, fieldName )
% % set vtk output
% if (obj.mesh.type == enumMeshDim.One)
%     obj.vtkOutput = VtkOutput1d(obj.casename, obj.Nfield, obj.timeInterval);
% elseif (obj.mesh.type == enumMeshDim.Two)
%     obj.vtkOutput = VtkOutput2d(obj.mesh, obj.casename, obj.Nfield, obj.timeInterval, varIndex);
% elseif (obj.mesh.type == enumMeshDim.Three)
%     obj.vtkOutput = VtkOutput3d(obj.mesh3d, obj.casename, obj.Nfield3d, obj.timeInterval3d);
% end

% obj.vtkOutput.initFromMesh( mesh );

obj.vtkOutput = VtkOutput( mesh, fieldname, casename, OutputFieldNum, dt, VarIndex );
obj.vtkOutput.initFromMesh( mesh );

% define dimension
dimTime = NdgNcDim('Nt', 0);
dimK = NdgNcDim('K', mesh.K);
dimNp = NdgNcDim('Np', mesh.cell.Np);
dimNfield = NdgNcDim('Nvar', obj.Nfield);
% dimNfield = NdgNcDim('Nvar', 6);
% define variable
varTime = NdgNcVar('time', dimTime, enumNcData.NC_DOUBLE );
varField = NdgNcVar('fphys', [dimNp, dimK, dimNfield, dimTime], enumNcData.NC_DOUBLE);


obj.ncid = zeros(size(filename));
obj.isOpen = false * ones(size(filename));
obj.fileOrder = 1;
obj.Numfile = numel(filename);
obj.filename = filename;
% define file
% obj.filename = [ obj.casename, '/', obj.casename, '.nc' ];
obj.ncfile = NdgNcFile( obj, ...
    [dimTime, dimK, dimNp, dimNfield], [varTime, varField]);

obj.StepPerFile = ceil(outputIntervalNum/numel(obj.filename));

% obj.ncfile.varIndex = varIndex;
% init file
for n = 1:numel(obj.filename)
    obj.defineIntoNetcdfFile(n);
end

% set properties
obj.timeVarableId = varTime.id;
obj.fieldVarableId = varField.id;
% obj.mesh = mesh;

end