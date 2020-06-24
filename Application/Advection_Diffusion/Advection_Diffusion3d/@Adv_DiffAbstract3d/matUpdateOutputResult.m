function matUpdateOutputResult( obj, time, fphys3d )
%MATUPDATEOUTPUTRESULT Summary of this function goes here
%   Detailed explanation goes here

for m = obj.Nmesh
    obj.outputFile3d(m).outputIntervalResult( time, fphys3d{m}(:,:,obj.outputFile3d(m).varIndex) );
end

end

