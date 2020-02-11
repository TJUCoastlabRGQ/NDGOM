function initPhysFromOptions( obj, mesh2d, mesh3d )

initPhysFromOptions@SWEAbstract3d( obj, mesh2d, mesh3d );

[ obj.advectionSolver, obj.viscositySolver, obj.numfluxSolver, obj.PCESolver2d ] = initSolver( obj );

obj.coriolisSolver = initcoriolisSolver3d(obj);

end