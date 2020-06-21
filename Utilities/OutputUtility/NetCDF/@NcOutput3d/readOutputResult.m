function [ field2d, field3d ] = readOutputResult( obj, step )


if ( step > obj.outputStep )
    error( ['The output step number is ',  num2str(obj.outputStep), ...
        ' less than input step number ', num2str(step) , '.\n'] )
end

Np = obj.mesh3d.cell.Np;
K = obj.mesh3d.K;
field3d = ncread( obj.ncfile.NcOutPut.fileName3d{1}, 'fphys3d', [1, 1, 1, step], [Np, K, obj.Nfield3d, 1]);
field2d = readOutputResult@NcOutput( obj, step );
end% func