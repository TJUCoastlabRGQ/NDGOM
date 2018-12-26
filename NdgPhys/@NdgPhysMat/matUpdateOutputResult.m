function matUpdateOutputResult( obj, time, fphys )

for m = 1:obj.Nmesh
%     varFieldIndex = [1 2 3 4 5 6];
    obj.outputFile(m).outputIntervalResult( time, fphys{m}(:,:,obj.varFieldIndex) );
%     obj.outputFile(m).outputIntervalResult( time, fphys{m}(:,:,varFieldIndex) );
end

end