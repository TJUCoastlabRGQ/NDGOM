function outputResult( obj, time, field )
    % output time
    startInd = obj.outputStep;
    countInd = 1;
    netcdf.putVar(obj.ncfile.ncid(obj.ncfile.fileOrder), obj.timeVarableId, startInd, countInd, time);
    
    % output physical field
    startInd = [ 0, 0, 0, obj.outputStep ];
    countInd = [ obj.mesh.cell.Np, obj.mesh.K, obj.Nfield, 1 ];
    netcdf.putVar(obj.ncfile.ncid(obj.ncfile.fileOrder), obj.fieldVarableId, startInd, countInd, field);
    
    if obj.outputStep == ( obj.ncfile.StepPerFile - 1 ) && obj.ncfile.fileOrder ~= obj.ncfile.Numfile
        obj.ncfile.closeNetcdfFile(obj.ncfile.fileOrder);
        obj.ncfile.fileOrder = obj.ncfile.fileOrder + 1;
        obj.outputStep = 0;
    else
    % increase output step num
    obj.outputStep = obj.outputStep + 1;        
    end
end
