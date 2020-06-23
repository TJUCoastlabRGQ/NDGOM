function writeOutputResultAtStepToVtk( obj, step )
%WRITEOUTPUTRESULTTOVTK Summary of this function goes here
%   Detailed explanation goes here
time = ncread(obj.fileName{1},'time');
fprintf('===============================\n');
fprintf('Time of the current output is:%f\n',time(step));
fprintf('Date has been written into the VTK file of order:%d\n',obj.vtkOutput.outputStep);
fprintf('===============================\n');
field = obj.readOutputResult( step );
obj.vtkOutput.outputResult( step, field );

end

