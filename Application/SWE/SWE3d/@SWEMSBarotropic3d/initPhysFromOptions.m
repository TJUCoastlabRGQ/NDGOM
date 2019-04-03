function initPhysFromOptions(obj, mesh2d, mesh3d)

initPhysFromOptions@SWEAbstract3d( obj, mesh2d, mesh3d );

[ obj.fphys2d, obj.fphys ] = obj.setInitialField;

obj.Solver2d.meshUnion = mesh2d;
obj.Solver2d.Nmesh = numel(mesh2d);
for m = 1:obj.Solver2d.Nmesh
    Np = obj.Solver2d.meshUnion(m).cell.Np;
    K = obj.Solver2d.meshUnion(m).K;
    obj.Solver2d.frhs{m} = zeros(Np, K, obj.Nvar2d);
    if ~isempty( mesh2d(m).BoundaryEdge )
        Nfp = mesh2d(m).BoundaryEdge.Nfp;
        Ne = mesh2d(m).BoundaryEdge.Ne;
        obj.Solver2d.fext{m} = zeros( Nfp, Ne, obj.Nfield );
    end
end

[ obj.Solver2d.advectionSolver, obj.Solver2d.viscositySolver ] = initSolver2d( obj );

if obj.option.isKey('CFL')
    obj.Solver2d.cfl = obj.getOption('CFL');
elseif obj.option.isKey('cfl')
    obj.Solver2d.cfl = obj.getOption('cfl');
elseif obj.option.isKey('Cfl')
    obj.Solver2d.cfl = obj.getOption('Cfl');
else
    obj.Solver2d.cfl = 1;
end

[ obj.Solver2d.limiter ] = initSlopeLimiter2d(obj);

obj.Solver2d.hmin = obj.hmin;
obj.Solver2d.matUpdateWetDryState( obj.fphys2d );

obj.Solver2d.zGrad = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    mesh = obj.mesh2d(m);
    
    obj.Solver2d.zGrad{m} = zeros( mesh2d.cell.Np, mesh2d.K, 2 );
    zr = mesh2d.cell.Dr * obj.fphys2d{m}(:,:,5);
    zs = mesh2d.cell.Ds * obj.fphys2d{m}(:,:,5);
    obj.Solver2d.zGrad{m}(:,:,1) = mesh.rx .* zr + mesh.sx .* zs;
    obj.Solver2d.zGrad{m}(:,:,2) = mesh.ry .* zr + mesh.sy .* zs;
end
end