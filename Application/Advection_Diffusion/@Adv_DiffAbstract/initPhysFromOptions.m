function initPhysFromOptions( obj, mesh )

% set mesh object
obj.meshUnion = mesh;
obj.Nmesh = numel( mesh );

% set the physical field for the NdgPhysMat solver
for m = 1:obj.Nmesh
    Np = obj.meshUnion(m).cell.Np;
    K = obj.meshUnion(m).K;
    obj.frhs{m} = zeros( Np, K, obj.Nvar );
    %> This is used to contain the right hand side corresponding to the implicit discretization of vertical diffusion term
    %> Only one mesh considered
    %     obj.ImplicitRHS3d = zeros( Np, K, obj.Nvar );
    if ~isempty( mesh(m).BoundaryEdge )
        Nfp = obj.meshUnion(m).BoundaryEdge.Nfp;
        Ne = obj.meshUnion(m).BoundaryEdge.Ne;
        obj.BoundaryEdgefext{m} = zeros( Nfp, Ne, obj.Nfield );
        %> The following two variables are used when calculating the horizontal diffusion term
        obj.BoundaryEdgefm{m} = zeros( Nfp, Ne, obj.Nfield );
        obj.BoundaryEdgefp{m} = zeros( Nfp, Ne, obj.Nfield );
    end
    
    %> These public variable are used to avoid the calculation of these terms repeatly during computation
    Nfp = obj.meshUnion(m).InnerEdge.Nfp;
    Ne = obj.meshUnion(m).InnerEdge.Ne;
    %> The following two variables are used when calculating the horizontal diffusion term
    obj.InnerEdgefm{m} = zeros( Nfp, Ne, obj.Nfield );
    obj.InnerEdgefp{m} = zeros( Nfp, Ne, obj.Nfield );
    
end

% set option
obj.option = obj.setOption( obj.option );

% [ obj.advectionSolver, obj.HorizontalEddyViscositySolver, obj.VerticalEddyViscositySolver ] = obj.initSolver();

[ obj.fphys ] = obj.setInitialField;


finalTime = obj.getOption('finalTime');
for m = 1:obj.Nmesh
    %     obj.fext{m} = obj.getExtFunc( obj.meshUnion(m), finalTime );
    obj.fext{m} = obj.getExtFunc( obj.meshUnion(m).x, obj.meshUnion(m).y, finalTime );
end
obj.outputFile = obj.matInitOutput(mesh, obj.fieldName);
end