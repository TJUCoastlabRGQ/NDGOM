function initFromMesh( obj, filename2d, filename3d, outputIntervalNum, varIndex2d, varIndex3d )

% set vtk output

obj.vtkOutput = VtkOutput2d(obj.mesh, obj.casename, obj.Nfield, obj.timeInterval, varIndex2d);
obj.vtkOutput3d = VtkOutput3d(obj.mesh3d, obj.casename, obj.Nfield3d, obj.timeInterval, varIndex3d);
% obj.vtkOutput.initFromMesh( mesh );

dimTime = NdgNcDim('Nt', 0);
dimK2 = NdgNcDim('K2d', obj.mesh.K);
dimNp2 = NdgNcDim('Np2d', obj.mesh.cell.Np);
dimNfield2 = NdgNcDim('Nfield2d', numel( varIndex2d ) );

dimK3 = NdgNcDim('K3d', obj.mesh3d.K );
dimNp3 = NdgNcDim('Np3d', obj.mesh3d.cell.Np );
dimNfield3 = NdgNcDim('Nfield3d', numel( varIndex3d ) );

% define variable
varTime = NdgNcVar('time', dimTime, enumNcData.NC_DOUBLE );
varField2 = NdgNcVar('fphys2d', ...
    [ dimNp2, dimK2, dimNfield2, dimTime], ...
    enumNcData.NC_DOUBLE);

varField3 = NdgNcVar('fphys3d', ...
    [ dimNp3, dimK3, dimNfield3, dimTime], ...
    enumNcData.NC_DOUBLE);

obj.ncfile = NdgNcFile( obj, ...
    [dimTime, dimK2, dimNp2, dimNfield2], ...
    [varTime, varField2]);

obj.ncfile3d = NdgNcFile( obj, ...
    [dimTime,dimK3, dimNp3, dimNfield3], ...
    [varTime, varField3]);
%

obj.ncid = zeros(size(filename2d));obj.ncid3d = zeros(size(filename3d));
obj.isOpen = false * ones(size(filename2d));obj.isOpen3d = false * ones(size(filename3d));
obj.fileOrder = 1; obj.fileOrder3d = 1;
obj.Numfile = numel(filename2d);
obj.fileName = filename2d;obj.fileName3d = filename3d;        

if floor(outputIntervalNum/numel(obj.fileName))<1
    error( 'Too many output nc file!' );
else
    obj.StepPerFile = floor(outputIntervalNum/outputIntervalNum);
end
%
% obj.ncfile.varIndex = varIndex;
% % init file
for n = 1:numel(obj.fileName)
    obj.defineIntoNetcdfFile(n);
end
% % set properties
obj.timeVarableId = varTime.id;
obj.timeVarableId3d = varTime.id;
obj.fieldVarableId = varField2.id;
obj.fieldVarableId3d = varField3.id;
end