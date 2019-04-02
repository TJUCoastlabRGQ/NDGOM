function matUpdateFinalResult( obj, time, fphys2d, fphys3d )
%MATUPDATEOUTPUTRESULT Summary of this function goes here
%   Detailed explanation goes here

for m = obj.Nmesh
    obj.outputFile(m).outputFinalResult( obj, time, fphys2d{m}, fphys3d{m} );
    obj.outputFile(m).ncfile.closeNetcdfFile(obj.outputFile(m).ncfile.fileOrder);
end

end