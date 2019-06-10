function initPhysFromOptions( obj, mesh )

% call the superclass methods
initPhysFromOptions@NdgPhysMat( obj, mesh );
obj.matUpdateWetDryState( obj.fphys );

obj.zGrad = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    mesh = obj.meshUnion(m);
    obj.zGrad{m} = zeros( mesh.cell.Np, mesh.K, 2 );
    zr = mesh.cell.Dr * obj.fphys{m}(:,:,3);
    obj.zGrad{m}(:,:,1) = mesh.rx .* zr;
end

%Friction Term
if obj.option.isKey('FrictionType') % the option exist
    switch obj.getOption('FrictionType')
        case SWEFrictionType.None
            obj.frictionSolver = NonFrictionTermSolver();
        case SWEFrictionType.Linear
            t = obj.getOption('FrictionCoefficient_r');
            obj.frictionSolver = LinearFrictionTermSolver1d(t);
        case SWEFrictionType.Manning
            n = obj.getOption('FrictionCoefficient_n');
            obj.frictionSolver = ManningFrictionSolver1d(n);
    end
else % the option does not exist
    obj.frictionSolver = NonFrictionTermSolver();
end

[ obj.numfluxSolver, obj.surfluxSolver, obj.volumefluxSolver ] = initFluxSolver( obj );

%Set NonhydrostaticSolver
[obj.NonhydrostaticSolver] = initNonhydrostaticSolver(obj);

if obj.option.isKey('SWELimiterType')
    switch obj.getOption('SWELimiterType')
        case SWELimiterType.OnDepth
            obj.limiterSolver = SWEDepthLimiter1d();
        case SWELimiterType.OnElevation
            obj.limiterSolver = SWEElevationLimiter1d();
    end
else
    obj.limiterSolver = SWEDepthLimiter1d();
end

end% func

function [NonhydrostaticSolver] = initNonhydrostaticSolver(obj)
if obj.option.isKey('nonhydrostaticType')
    switch obj.getOption('nonhydrostaticType')
        case enumNonhydrostaticType.Hydrostatic
            NonhydrostaticSolver = NdghydrostaticSolver1d(obj);
        case enumNonhydrostaticType.Nonhydrostatic
            % if Nonhydrostatic Solver included , the output variable is [ h, hu, hv, hw ]
            switch obj.getOption('integralType')
                case enumDiscreteIntegral.GaussQuadrature
                    NonhydrostaticSolver = NdgGaussQuadNonhydrostaticSolver1d(obj, obj.meshUnion);
                case enumDiscreteIntegral.QuadratureFree
                    NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver1d(obj);
            end
    end% switch
else
    NonhydrostaticSolver = NdghydrostaticSolver1d( obj );
end
end

