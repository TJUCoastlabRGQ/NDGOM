function writeOutputResultAtTimePointToVtk( obj, timepoint )
%WRITEOUTPUTRESULTTOVTK Summary of this function goes here
%   Detailed explanation goes here
time = ncread(obj.fileName3d{1},'time');
[~,Index] = sort(abs(timepoint-time));		   
fprintf('===============================\n');
fprintf('Time of the current output is:%f\n',time(Index(1)));
fprintf('Date has been written into the VTK file of order:%d\n',obj.vtkOutput3d.outputStep3d);
fprintf('===============================\n');
[ field2d, field3d ] = obj.readOutputResult( Index(1) );
obj.vtkOutput3d.outputStepResult( field2d, field3d );
end

