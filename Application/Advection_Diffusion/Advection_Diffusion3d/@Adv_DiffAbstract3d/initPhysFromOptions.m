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
        obj.BoundaryEdgefext3d{m} = zeros( Nfp, Ne, obj.Nfield );
        %> The following two variables are used when calculating the horizontal diffusion term
        obj.BoundaryEdgefm3d{m} = zeros( Nfp, Ne, obj.Nvar );
        obj.BoundaryEdgefp3d{m} = zeros( Nfp, Ne, obj.Nvar );
    end
    if ~isempty( mesh3d(m).SurfaceBoundaryEdge )
        Nfp = obj.meshUnion(m).SurfaceBoundaryEdge.Nfp;
        Ne = obj.meshUnion(m).SurfaceBoundaryEdge.Ne;
        obj.SurfaceBoundaryfext3d{m} = zeros( Nfp, Ne, obj.Nfield );
        %> The following two variables are used when calculating the horizontal diffusion term
        obj.SurfaceBoundaryEdgefm3d{m} = zeros( Nfp, Ne, obj.Nvar );
        obj.SurfaceBoundaryEdgefp3d{m} = zeros( Nfp, Ne, obj.Nvar );
    end
    if ~isempty( mesh3d(m).BottomBoundaryEdge )
        Nfp = obj.meshUnion(m).BottomBoundaryEdge.Nfp;
        Ne = obj.meshUnion(m).BottomBoundaryEdge.Ne;
        obj.BottomBoundaryEdgefext3d{m} = zeros( Nfp, Ne, obj.Nfield );
        %> The following two variables are used when calculating the horizontal diffusion term
        obj.BottomBoundaryEdgefm3d{m} = zeros( Nfp, Ne, obj.Nvar );
        obj.BottomBoundaryEdgefp3d{m} = zeros( Nfp, Ne, obj.Nvar );
    end    
    
    %> These public variable are used to avoid the calculation of these terms repeatly during computation
    Nfp = obj.meshUnion(m).InnerEdge.Nfp;
    Ne = obj.meshUnion(m).InnerEdge.Ne;
    %> The following two variables are used when calculating the horizontal diffusion term
    obj.InnerEdgefm3d{m} = zeros( Nfp, Ne, obj.Nfield );
    obj.InnerEdgefp3d{m} = zeros( Nfp, Ne, obj.Nfield );
   
end

% set option
obj.option = obj.setOption( obj.option );

[ obj.advectionSolver, obj.HorizontalEddyViscositySolver, obj.VerticalEddyViscositySolver ] = obj.initSolver();

%> set the final time
obj.ftime = obj.getOption('finalTime');

[ obj.fphys ] = obj.setInitialField;

% we currently only consider one mesh only
obj.SurfBoundNewmannDate = zeros(mesh2d(1).cell.Np, mesh2d(1).K, obj.Nvar);
obj.BotBoundNewmannDate = zeros(mesh2d(1).cell.Np, mesh2d(1).K, obj.Nvar);

obj.outputFile3d = obj.matInitOutput(mesh3d, obj.fieldName3d);
end