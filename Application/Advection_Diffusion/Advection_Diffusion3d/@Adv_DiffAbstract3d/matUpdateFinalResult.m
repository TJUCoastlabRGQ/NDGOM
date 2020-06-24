function matUpdateFinalResult( obj, time, fphys3d )
%MATUPDATEOUTPUTRESULT Summary of this function goes here
%   Detailed explanation goes here

for m = obj.Nmesh
    obj.outputFile3d(m).outputFinalResult( time, fphys3d{m}(:,:,obj.outputFile3d(m).varIndex) );
%     obj.outputFile(m).closeNetcdfFile(obj.outputFile(m).fileOrder);
end

end