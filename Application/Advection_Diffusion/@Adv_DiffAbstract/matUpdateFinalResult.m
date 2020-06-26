function matUpdateFinalResult( obj, time, fphys )
%MATUPDATEOUTPUTRESULT Summary of this function goes here
%   Detailed explanation goes here

for m = obj.Nmesh
    obj.outputFile(m).outputFinalResult( time, fphys{m}(:,:,obj.outputFile(m).varIndex) );
%     obj.outputFile(m).closeNetcdfFile(obj.outputFile(m).fileOrder);
end

end