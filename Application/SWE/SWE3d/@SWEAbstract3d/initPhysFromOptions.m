function initPhysFromOptions( obj, mesh2d, mesh3d )

% set mesh object
obj.mesh2d = mesh2d;
obj.meshUnion = mesh3d;
obj.Nmesh = numel( mesh3d );

% set the physical field for the NdgPhysMat solver
for m = 1:obj.Nmesh
    Np = obj.meshUnion(m).cell.Np;
    K = obj.meshUnion(m).K;
    obj.frhs{m} = zeros( Np, K, obj.Nvar );
    %> This is used to contain the right hand side corresponding to the implicit discretization of vertical diffusion term
    %> Only one mesh considered
%     obj.ImplicitRHS3d = zeros( Np, K, obj.Nvar );
    if ~isempty( mesh3d(m).BoundaryEdge )
        Nfp = obj.meshUnion(m).BoundaryEdge.Nfp;
        Ne = obj.meshUnion(m).BoundaryEdge.Ne;
        obj.fext3d{m} = zeros( Nfp, Ne, obj.Nfield );
        %> The following variable are used when calculating the horizontal diffusion term
        obj.BoundaryEdgefp{m} = zeros( Nfp, Ne, obj.Nvar );
    end
    %> These public variable are used to avoid the calculation of these terms repeatly during computation
    
    Np = obj.mesh2d(m).cell.Np;
    K = obj.mesh2d(m).K;
    obj.frhs2d{m} = zeros( Np, K, obj.Nvar2d );
    if ~isempty( mesh2d(m).BoundaryEdge )
        Nfp = obj.mesh2d(m).BoundaryEdge.Nfp;
        Ne = obj.mesh2d(m).BoundaryEdge.Ne;
        obj.fext2d{m} = zeros( Nfp, Ne, obj.Nfield2d );
    end
end

% set option
obj.option = obj.setOption( obj.option );

%> set the final time
obj.ftime = obj.getOption('finalTime');
if obj.option.isKey('CFL')
    obj.cfl = obj.getOption('CFL');
elseif obj.option.isKey('cfl')
    obj.cfl = obj.getOption('cfl');
elseif obj.option.isKey('Cfl')
    obj.cfl = obj.getOption('Cfl');
else
    obj.cfl = 1;
end

%> bottom friction coefficient
obj.Cf =  cell(obj.Nmesh);

[ obj.fphys2d, obj.fphys ] = obj.setInitialField;
obj.matUpdateWetDryState(obj.fphys2d);
for n = 1:obj.Nmesh
    obj.fphys2d{n}(:,:,5) = ( obj.mesh2d.rx .* ( obj.mesh2d.cell.Dr * obj.fphys2d{n}(:,:,4) ) ) + ...
        ( obj.mesh2d.sx .* ( obj.mesh2d.cell.Ds * obj.fphys2d{n}(:,:,4) ) );
    obj.fphys2d{n}(:,:,6) = ( obj.mesh2d.ry .* ( obj.mesh2d.cell.Dr * obj.fphys2d{n}(:,:,4) ) ) + ...
        ( obj.mesh2d.sy .* ( obj.mesh2d.cell.Ds * obj.fphys2d{n}(:,:,4) ) );
    obj.fphys2d{n}(:, :, 2) = obj.meshUnion(n).VerticalColumnIntegralField( obj.fphys{n}(:, :, 1) );
    obj.fphys2d{n}(:, :, 3) = obj.meshUnion(n).VerticalColumnIntegralField( obj.fphys{n}(:, :, 2) );
    %> H in extended three dimensional fields
    obj.fphys{n}(:, :, 4) = obj.meshUnion(n).Extend2dField( obj.fphys2d{n}(:, :, 1) );
    %> the vertical viscosity term is initialized to be 10^(-2)
    obj.fphys{n}(:,:,5) = 0.0001* ones(size(obj.meshUnion(n).x));
    %> Initially, the bottom roughness is not considered
    obj.Cf{n} =  0.0025/1000 * ones(size(obj.mesh2d(n).y));
    %> Z in extended three dimensional fields
    obj.fphys{n}(:,:,6) = obj.meshUnion(n).Extend2dField( obj.fphys2d{n}(:, :, 4) );
    %> '$\eta$' in extended three dimensional fields
    obj.fphys{n}(:,:,7) = obj.fphys{n}(:,:,4) + obj.fphys{n}(:,:,6);
    %> Zx in extended three dimensional fields
    obj.fphys{n}(:,:,8) = obj.meshUnion(n).Extend2dField( obj.fphys2d{n}(:, :, 5) );
    %> Zy in extended three dimensional fields
    obj.fphys{n}(:,:,9) = obj.meshUnion(n).Extend2dField( obj.fphys2d{n}(:, :, 6) );
end

% we currently only consider one mesh only
obj.SurfBoundNewmannDate = zeros(mesh2d(1).cell.Np, mesh2d(1).K, obj.Nvar);
obj.BotBoundNewmannDate = zeros(mesh2d(1).cell.Np, mesh2d(1).K, obj.Nvar);

obj.matInitEddyViscositySolver( );
% Setup the output NetCDF file object
% initOutput( obj, mesh2d, mesh3d );
obj.outputFile3d = obj.matInitOutput(mesh3d, obj.fieldName3d);
obj.outputFile2d = obj.matInitOutput(mesh2d, obj.fieldName2d);
end