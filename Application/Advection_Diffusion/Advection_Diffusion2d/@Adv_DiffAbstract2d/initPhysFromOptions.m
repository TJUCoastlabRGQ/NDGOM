function initPhysFromOptions( obj, mesh )

initPhysFromOptions@Adv_DiffAbstract( obj, mesh );

[ obj.advectionSolver, obj.HorizontalEddyViscositySolver, obj.VerticalEddyViscositySolver ] = obj.initSolver();

end