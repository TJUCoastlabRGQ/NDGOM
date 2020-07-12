function writeOutputResultToVtk( obj )
%WRITEOUTPUTRESULTTOVTK Summary of this function goes here
%   Detailed explanation goes here
time = ncread(obj.filename{1},'time');
obj.OutputTime = time;
fprintf('Time of the output is:\n');
disp( obj.OutputTime );
fprintf('===============================\n');
for step = 1:numel(time)
    field = obj.readOutputResult( step );
    obj.vtkOutput.outputResult( field );
end

end