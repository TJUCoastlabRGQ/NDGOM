function initPhysFromOptions( obj, mesh )

initPhysFromOptions@Adv_DiffAbstract( obj, mesh );

[ obj.advectionSolver, obj.HorizontalEddyViscositySolver, obj.VerticalEddyViscositySolver ] = obj.initSolver();

obj.fext2d = cell(obj.Nmesh);
end