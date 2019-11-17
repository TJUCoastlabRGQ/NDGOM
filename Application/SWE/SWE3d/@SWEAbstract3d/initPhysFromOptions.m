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

% set option
obj.option = obj.setOption( obj.option );

obj.EddyViscositySolver = obj.matInitEddyViscositySolver( );

%> wind stress term
obj.WindTaux =  cell(obj.Nmesh);
obj.WindTauy =  cell(obj.Nmesh);
%> bottom friction coefficient
obj.Cf =  cell(obj.Nmesh);

[ obj.fphys2d, obj.fphys ] = obj.setInitialField;

for n = 1:obj.Nmesh
    obj.fphys2d{n}(:,:,5) = ( mesh2d.rx .* ( mesh2d.cell.Dr * obj.fphys2d{n}(:,:,4) ) ) + ...
        ( mesh2d.sx .* ( mesh2d.cell.Ds * obj.fphys2d{n}(:,:,4) ) );
    obj.fphys2d{n}(:,:,6) = ( mesh2d.ry .* ( mesh2d.cell.Dr * obj.fphys2d{n}(:,:,4) ) ) + ...
        ( mesh2d.sy .* ( mesh2d.cell.Ds * obj.fphys2d{n}(:,:,4) ) );
    %> the wind stress is initialized to be zero
    obj.WindTaux{n} = zeros(size(mesh2d(n).x));
    obj.WindTauy{n} = zeros(size(mesh2d(n).y));
    %> H in extended three dimensional fields
    obj.fphys{n}(:, :, 4) = mesh3d(n).Extend2dField( obj.fphys2d{n}(:, :, 1) );
    %> the vertical viscosity term is initialized to be 10^(-4)
    obj.fphys{n}(:,:,5) = 1e-4* ones(size(mesh3d(n).x)); 
    %> Initially, the bottom roughness is not considered
    obj.Cf{n} =  0.0025/1000 * ones(size(mesh2d(n).y)); 
    %> Z in extended three dimensional fields
    obj.fphys{n}(:,:,6) = mesh3d(n).Extend2dField( obj.fphys2d{n}(:, :, 4) );
    %> '$\eta$' in extended three dimensional fields
    obj.fphys{n}(:,:,7) = obj.fphys{n}(:,:,4) + obj.fphys{n}(:,:,6);
    %> Zx in extended three dimensional fields
    obj.fphys{n}(:,:,8) = mesh3d(n).Extend2dField( obj.fphys2d{n}(:, :, 5) );
    %> Zy in extended three dimensional fields
    obj.fphys{n}(:,:,9) = mesh3d(n).Extend2dField( obj.fphys2d{n}(:, :, 6) );
end
% Setup the output NetCDF file object
initOutput( obj, mesh2d, mesh3d );
end