function writeOutputResultAtTimePointToVtk( obj, timepoint )
%WRITEOUTPUTRESULTTOVTK Summary of this function goes here
%   Detailed explanation goes here
time = ncread(obj.fileName{1},'time');
[~,Index] = sort(abs(timepoint-time));		   
fprintf('===============================\n');
fprintf('Time of the current output is:%f\n',time(Index(1)));
fprintf('Date has been written into the VTK file of order:%d\n',obj.vtkOutput.outputStep);
fprintf('===============================\n');
field = obj.readOutputResult( Index(1) );
obj.vtkOutput3d.outputStepResult( field );

end

