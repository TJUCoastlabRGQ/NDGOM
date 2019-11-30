function outputResult( obj, time, field )
    % output time
    startInd = obj.outputStep;
    countInd = 1;
    netcdf.putVar(obj.ncid(obj.fileOrder), obj.timeVarableId, startInd, countInd, time);
    
    % output physical field
    startInd = [ 0, 0, 0, obj.outputStep ];
    countInd = [ obj.mesh.cell.Np, obj.mesh.K, obj.Nfield, 1 ];
    netcdf.putVar(obj.ncid(obj.fileOrder), obj.fieldVarableId, startInd, countInd, field);
    
    if obj.outputStep == ( obj.StepPerFile - 1 ) && obj.fileOrder ~= obj.Numfile
        obj.isOpen(obj.fileOrder) = false;
        netcdf.close(obj.ncid(obj.fileOrder));
        obj.fileOrder = obj.fileOrder + 1;
        obj.outputStep = 0;
    else
    % increase output step num
    obj.outputStep = obj.outputStep + 1;        
    end
end
