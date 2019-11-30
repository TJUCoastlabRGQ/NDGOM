function matUpdateOutputResult( obj, time, fphys )

for m = 1:obj.Nmesh
    obj.outputFile(m).outputIntervalResult( time, fphys{m}(:,:,obj.outputFile(m).varIndex) );
end

end