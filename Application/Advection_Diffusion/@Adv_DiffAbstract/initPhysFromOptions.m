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
    
end

% set option
obj.option = obj.setOption( obj.option );

% [ obj.advectionSolver, obj.HorizontalEddyViscositySolver, obj.VerticalEddyViscositySolver ] = obj.initSolver();

[ obj.fphys ] = obj.setInitialField;


finalTime = obj.getOption('finalTime');
for m = 1:obj.Nmesh
        obj.ExactValue{m} = obj.getExtFunc( obj.meshUnion(m), finalTime );
end
obj.outputFile = obj.matInitOutput(mesh, obj.fieldName);
end