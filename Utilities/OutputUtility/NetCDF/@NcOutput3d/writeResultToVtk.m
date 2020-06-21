function writeResultToVtk(obj, step, field)
%myFun - Description
%
% Syntax: output = myFun(input)
%
% Long description

writeResultToVtk@NcOutput(obj, step, field);
obj.vtkOutput.outputStepResult( step, field );
end