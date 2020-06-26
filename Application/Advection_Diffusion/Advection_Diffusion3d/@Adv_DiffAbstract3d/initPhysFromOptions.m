function initPhysFromOptions( obj, mesh2d, mesh3d )

initPhysFromOptions@Adv_DiffAbstract( obj, mesh3d);
% set mesh object
obj.mesh2d = mesh2d;

% set the physical field for the NdgPhysMat solver
for m = 1:obj.Nmesh
    if ~isempty( mesh3d(m).SurfaceBoundaryEdge )
        Nfp = obj.meshUnion(m).SurfaceBoundaryEdge.Nfp;
        Ne = obj.meshUnion(m).SurfaceBoundaryEdge.Ne;
        obj.SurfaceBoundaryfext{m} = zeros( Nfp, Ne, obj.Nfield );
        %> The following two variables are used when calculating the horizontal diffusion term
        obj.SurfaceBoundaryEdgefm{m} = zeros( Nfp, Ne, obj.Nvar );
        obj.SurfaceBoundaryEdgefp{m} = zeros( Nfp, Ne, obj.Nvar );
    end
    if ~isempty( mesh3d(m).BottomBoundaryEdge )
        Nfp = obj.meshUnion(m).BottomBoundaryEdge.Nfp;
        Ne = obj.meshUnion(m).BottomBoundaryEdge.Ne;
        obj.BottomBoundaryEdgefext{m} = zeros( Nfp, Ne, obj.Nfield );
        %> The following two variables are used when calculating the horizontal diffusion term
        obj.BottomBoundaryEdgefm{m} = zeros( Nfp, Ne, obj.Nvar );
        obj.BottomBoundaryEdgefp{m} = zeros( Nfp, Ne, obj.Nvar );
    end 
end

[ obj.advectionSolver, obj.HorizontalEddyViscositySolver, obj.VerticalEddyViscositySolver ] = obj.initSolver();

% we currently only consider one mesh only
obj.SurfBoundNewmannDate = zeros(mesh2d(1).cell.Np, mesh2d(1).K, obj.Nvar);
obj.BotBoundNewmannDate = zeros(mesh2d(1).cell.Np, mesh2d(1).K, obj.Nvar);


obj.fext3d = cell(obj.Nmesh);
end