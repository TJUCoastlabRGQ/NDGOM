function outputResult( obj, time, field2d, field3d )

outputResult@NcOutput(obj, time, field2d);

% output physical field
startInd = [ 0, 0, 0, obj.outputStep3d ];
countInd = [ obj.mesh3d.cell.Np, obj.mesh3d.K, obj.Nfield3d, 1 ];
netcdf.putVar(obj.ncid3d(obj.fileOrder3d), obj.fieldVarableId3d, startInd, countInd, field3d);

if obj.outputStep3d == ( obj.StepPerFile - 1 ) && obj.fileOrder3d ~= obj.Numfile
    obj.isOpen3d(obj.fileOrder) = false;
    netcdf.close(obj.ncid3d(obj.fileOrder3d));
    obj.fileOrder3d = obj.fileOrder3d + 1;
    obj.outputStep3d = 0;
else
    % increase output step num
    obj.outputStep3d = obj.outputStep3d + 1;
end
end