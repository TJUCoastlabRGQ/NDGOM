function initPhysFromOptions( obj, mesh2d, mesh3d )

initPhysFromOptions@SWEAbstract3d( obj, mesh2d, mesh3d );

[ obj.advectionSolver, obj.viscositySolver, obj.windSolver, obj.numfluxSolver, obj.frictionSolver, obj.PCESolver2d ] = initSolver( obj );

end