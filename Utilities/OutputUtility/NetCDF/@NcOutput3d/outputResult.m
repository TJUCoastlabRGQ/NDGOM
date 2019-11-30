function outputResult( obj, time, field2d, field3d )
    % output time
%     outputObj = matInitOutput@NdgPhysMat(obj);
    outputResult@NcOutput(obj, time, field2d);
    startInd = obj.outputStep3d;
    countInd = 1;
    netcdf.putVar(obj.ncfile3d.ncid(obj.ncfile.fileOrder), obj.timeVarableId, startInd, countInd, time);
    
    % output physical field
    startInd = [ 0, 0, 0, obj.outputStep3d ];
    countInd = [ obj.mesh3d.cell.Np, obj.mesh3d.K, obj.Nfield3d, 1 ];
    netcdf.putVar(obj.ncfile3d.ncid(obj.ncfile3d.fileOrder), obj.fieldVarableId3d, startInd, countInd, field3d);
    
    if obj.outputStep == ( obj.ncfile3d.StepPerFile - 1 ) && obj.ncfile3d.fileOrder ~= obj.ncfile3d.Numfile
        obj.ncfile3d.closeNetcdfFile(obj.ncfile3d.fileOrder);
        obj.ncfile3d.fileOrder = obj.ncfile3d.fileOrder + 1;
        obj.outputStep = 0;
    else
    % increase output step num
    obj.outputStep3d = obj.outputStep3d + 1;        
    end
end