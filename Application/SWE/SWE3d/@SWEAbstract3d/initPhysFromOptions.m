function initPhysFromOptions( obj, mesh2d, mesh3d )

% set mesh object
obj.mesh2d = mesh2d;
obj.meshUnion = mesh3d;
obj.Nmesh = numel( mesh3d );

% set the physical field for the NdgPhysMat solver
for m = 1:obj.Nmesh
    Np = obj.mesh3d(m).cell.Np;
    K = obj.mesh3d(m).K;
    obj.frhs{m} = zeros( Np, K, obj.Nvar );
    if ~isempty( mesh3d(m).BoundaryEdge )
        Nfp = mesh3d(m).BoundaryEdge.Nfp;
        Ne = mesh3d(m).BoundaryEdge.Ne;
        obj.fext3d{m} = zeros( Nfp, Ne, obj.Nfield );
    end
    
    Np = obj.mesh2d(m).cell.Np;
    K = obj.mesh2d(m).K;
    obj.frhs2d{m} = zeros( Np, K, obj.Nvar2d );
    if ~isempty( mesh2d(m).BoundaryEdge )
        Nfp = mesh2d(m).BoundaryEdge.Nfp;
        Ne = mesh2d(m).BoundaryEdge.Ne;
        obj.fext2d{m} = zeros( Nfp, Ne, obj.Nfield2d );
    end
end
%> wind stress term
obj.Taux =  cell(obj.Nmesh);
obj.Tauy =  cell(obj.Nmesh);
%> bottom friction coefficient
obj.Cf =  cell(obj.Nmesh);
for n = 1:obj.Nmesh
    obj.Taux{n} = zeros(size(mesh2d(n).x));
    obj.Tauy{n} = zeros(size(mesh2d(n).y));
    obj.Cf{n} =  zeros(size(mesh2d(n).y));
end
% Setup the output NetCDF file object
initOutput( obj, mesh2d, mesh3d );
end