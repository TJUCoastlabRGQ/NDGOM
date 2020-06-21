function outputResult(obj, time, field2d, field3d)
%outputResult - Description
%
% Syntax: outputResult(obj, time, field)
%
% Long description
outputResult@VtkOutput2d(obj, time, field2d);    

filename = [ obj.casename, '/3d/', obj.casename, '.', ...
    num2str(obj.outputStep3d, '%04d'), '.vtk' ];
fp = fopen(filename, 'w');

fprintf(fp, '# vtk DataFile Version 2.0\n');
fprintf(fp, ['NDG-FEM ', obj.casename, ' Simulation\n']);
fprintf(fp, 'ASCII\n');
fprintf(fp, 'DATASET UNSTRUCTURED_GRID\n');

fprintf(fp, '\nPOINTS %d double\n', obj.Npoint3d);
fprintf(fp, '%12.20f  %12.20f  %12.20f\n', obj.Points3d);

fprintf(fp, '\nCELLS %d %d\n', obj.Ncell3d, obj.Ncon3d);
dataFormat = [num2str(obj.Np3d, '%d'), ' ', repmat('%d ', 1, obj.Np3d), '\n'];
fprintf(fp, dataFormat, obj.CellVertList3d);

fprintf(fp, '\nCELL_TYPES %d\n', obj.Ncell3d);
dataFormat = [repmat('%d ', 1, 12), '\n'];
fprintf(fp, dataFormat, obj.ctype3d);

fprintf(fp, '\n\nPOINT_DATA %d\n', obj.Npoint3d);

for i = 1:size(field3d,3)
    fprintf(fp, ['SCALARS  ', obj.varName3d{i}, '  double\n']);
    fprintf(fp, 'LOOKUP_TABLE default\n');
    fprintf(fp, '%20.12f \n', field3d(:,:,i));
end

fclose(fp);

obj.outputStep3d = obj.outputStep3d + 1;
end