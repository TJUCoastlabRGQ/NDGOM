function initPhysFromOptions( obj, mesh2d, mesh3d )

initPhysFromOptions@Adv_DiffAbstract( obj, mesh3d);
% set mesh object
obj.mesh2d = mesh2d;

[ obj.advectionSolver, obj.HorizontalEddyViscositySolver, obj.VerticalEddyViscositySolver ] = obj.initSolver();

% we currently only consider one mesh only
obj.SurfBoundNewmannDate = zeros(mesh2d(1).cell.Np, mesh2d(1).K, obj.Nvar);
obj.BotBoundNewmannDate = zeros(mesh2d(1).cell.Np, mesh2d(1).K, obj.Nvar);

    for m = 1:obj.Nmesh
        if ~isempty( mesh3d(m).SurfaceBoundaryEdge )
            Nfp = mesh3d(m).SurfaceBoundaryEdge.Nfp;
            Ne = mesh3d(m).SurfaceBoundaryEdge.Ne;
            obj.SurfEdgefext{m} = zeros( Nfp, Ne, obj.Nfield );
        end
        if ~isempty( mesh3d(m).BottomBoundaryEdge )
            Nfp = mesh3d(m).BottomBoundaryEdge.Nfp;
            Ne = mesh3d(m).BottomBoundaryEdge.Ne;
            obj.BotEdgefext{m} = zeros( Nfp, Ne, obj.Nfield );
        end        
    end
end