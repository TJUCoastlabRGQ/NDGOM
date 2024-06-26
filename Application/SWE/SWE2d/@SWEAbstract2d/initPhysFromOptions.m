function initPhysFromOptions( obj, mesh )

% call the superclass methods
initPhysFromOptions@NdgPhysMat( obj, mesh );
obj.matUpdateWetDryState( obj.fphys );

obj.zGrad = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    
    obj.zGrad{m} = zeros( mesh.cell.Np, mesh.K, 2 );
    zr = mesh.cell.Dr * obj.fphys{m}(:,:,4);
    zs = mesh.cell.Ds * obj.fphys{m}(:,:,4);
    obj.zGrad{m}(:,:,1) = mesh.rx .* zr + mesh.sx .* zs;
    obj.zGrad{m}(:,:,2) = mesh.ry .* zr + mesh.sy .* zs;
end

%Set NonhydrostaticSolver
[obj.NonhydrostaticSolver] = initNonhydrostaticSolver(obj);
%Coriolis Term
[ obj.coriolisSolver ] = initCoriolisSolver( obj );
%Friction Term
[ obj.frictionSolver ] = initFrictionSolver( obj );
%Wind Term
[ obj.windSolver ] = initWindSolver( obj );
[ obj.numfluxSolver, obj.surfluxSolver, obj.volumefluxSolver ] = initFluxSolver( obj );
[ obj.limiterSolver ] = initLimiter( obj );

    % Setup the output NetCDF file object
obj.outputFile = obj.matInitOutput(mesh, obj.fieldName);

obj.VerticalEddyViscositySolver = NdgNoneVertDiffSolver( obj );

end% func

function [NonhydrostaticSolver] = initNonhydrostaticSolver(obj)
if obj.option.isKey('nonhydrostaticType')
    switch obj.getOption('nonhydrostaticType')
        case enumNonhydrostaticType.Hydrostatic
            NonhydrostaticSolver = NdghydrostaticSolver2d(obj);
        case enumNonhydrostaticType.Nonhydrostatic
            % if Nonhydrostatic Solver included , the output variable is [ h, hu, hv, hw ]
            switch obj.getOption('integralType')
                case enumDiscreteIntegral.GaussQuadrature
                    NonhydrostaticSolver = NdgGaussQuadNonhydrostaticSolver2d(obj, obj.meshUnion);
                case enumDiscreteIntegral.QuadratureFree
                    NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver2d(obj);
            end
    end% switch
else
    NonhydrostaticSolver = NdghydrostaticSolver2d( obj );
end
end
