function initFromMesh( obj, filename2d, filename3d, outputIntervalNum, varIndex2d, varIndex3d )

% set vtk output

obj.vtkOutput = VtkOutput2d(obj.mesh, obj.casename, obj.Nfield, obj.timeInterval);
obj.vtkOutput3d = VtkOutput3d(obj.mesh3d, obj.casename, obj.Nfield3d, obj.timeInterval3d);
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

obj.ncfile = NdgNcFile( obj, filename2d, ...
    [dimTime, dimK2, dimNp2, dimNfield2], ...
    [varTime, varField2]);

obj.ncfile3d = NdgNcFile( obj, filename3d, ...
    [dimTime,dimK3, dimNp3, dimNfield3], ...
    [varTime, varField3]);
%
if floor(outputIntervalNum/numel(obj.ncfile.fileName))<1
    error( 'Too many output nc file!' );
else
    obj.ncfile.StepPerFile = floor(outputIntervalNum/outputIntervalNum);
    obj.ncfile3d.StepPerFile = floor(outputIntervalNum/outputIntervalNum);
end
%
% obj.ncfile.varIndex = varIndex;
% % init file
for n = 1:numel(obj.ncfile.fileName)
    obj.ncfile.defineIntoNetcdfFile(n);
    obj.ncfile3d.defineIntoNetcdfFile(n);
end
% % set properties
obj.timeVarableId = varTime.id;
% obj.timeVarableId3d = varTime.id;
obj.fieldVarableId = varField2.id;
obj.fieldVarableId3d = varField3.id;
end