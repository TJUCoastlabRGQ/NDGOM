function matUpdateOutputResult( obj, time, fphys2d, fphys3d )
%MATUPDATEOUTPUTRESULT Summary of this function goes here
%   Detailed explanation goes here

for m = obj.Nmesh
    obj.outputFile(m).outputResult( obj, time, fphys2d{m}, fphys3d{m} )
end

end
