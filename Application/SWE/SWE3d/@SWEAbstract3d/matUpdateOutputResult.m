function matUpdateOutputResult( obj, time, fphys2d, fphys3d )
%MATUPDATEOUTPUTRESULT Summary of this function goes here
%   Detailed explanation goes here

for m = obj.Nmesh
    obj.outputFile(m).outputIntervalResult( time, fphys2d{m}(:,:,obj.outputFile(m).varIndex), fphys3d{m}(:,:,obj.outputFile(m).varIndex3d) )
end

end

