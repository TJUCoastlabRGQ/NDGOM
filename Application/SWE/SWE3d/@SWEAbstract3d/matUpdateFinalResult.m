function matUpdateFinalResult( obj, time, fphys2d, fphys3d )
%MATUPDATEOUTPUTRESULT Summary of this function goes here
%   Detailed explanation goes here

for m = obj.Nmesh
    obj.outputFile2d(m).outputFinalResult( time, fphys2d{m}(:,:,obj.outputFile2d(m).varIndex) );
    obj.outputFile3d(m).outputFinalResult( time, fphys3d{m}(:,:,obj.outputFile3d(m).varIndex) );
%     obj.outputFile(m).closeNetcdfFile(obj.outputFile(m).fileOrder);
end

end