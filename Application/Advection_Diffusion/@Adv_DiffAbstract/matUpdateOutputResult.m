function matUpdateOutputResult( obj, time, fphys )
%MATUPDATEOUTPUTRESULT Summary of this function goes here
%   Detailed explanation goes here

for m = obj.Nmesh
    obj.outputFile(m).outputIntervalResult( time, fphys{m}(:,:,obj.outputFile(m).varIndex) );
end

end

