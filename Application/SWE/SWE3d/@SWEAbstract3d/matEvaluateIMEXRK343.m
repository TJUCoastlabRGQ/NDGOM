function matEvaluateIMEXRK343( obj )
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
time = obj.startTime;
ftime = obj.finalTime;

fphys2d = obj.fphys2d;
fphys = obj.fphys;
%> allocate space for the rhs to be stored
ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,4);
ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 4*obj.Nvar);
ImplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
SystemRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 2);
DiffusionCoefficient = fphys{1}(:,:,5)./fphys{1}(:,:,4).^2;
visual = Visual2d( obj.mesh2d );
dt = obj.dt;
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    %     dt = obj.matUpdateTimeInterval( fphys2d );
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,1);
    Tempfphys = fphys{1}(:,:,1:2);

    for intRK = 1:3
        tloc = time + c( intRK ) * dt;
        %>Actually, boundary condition need to be imposed here
        %         obj.matUpdateExternalField( tloc, fphys2d, fphys3d );
        %This part need to consider the impact of the fext3d, as this is needed when impose the three-dimensional boundary
        [ExplicitRHS2d(:,:,intRK), ExplicitRHS3d(:,:,intRK), ExplicitRHS3d(:,:,intRK+4)] = ...
            matCalculateExplicitRHSTerm(obj, fphys2d, fphys, obj.fext2d);
        
        SystemRHS(:,:,1) = Tempfphys(:,:,1) + dt * EXa(intRK+1,1)*ExplicitRHS3d(:,:,1)+ dt * EXa(intRK+1,2)*ExplicitRHS3d(:,:,2)+...
            dt * EXa(intRK+1,3)*ExplicitRHS3d(:,:,3) + dt * IMa(intRK,1)*ImplicitRHS3d(:,:,1) + dt * IMa(intRK,2)*ImplicitRHS3d(:,:,2)+...
            dt * IMa(intRK,3)*ImplicitRHS3d(:,:,3);
        
        SystemRHS(:,:,2) = Tempfphys(:,:,2) + dt * EXa(intRK+1,1)*ExplicitRHS3d(:,:,5)+ dt * EXa(intRK+1,2)*ExplicitRHS3d(:,:,6)+...
            dt * EXa(intRK+1,3)*ExplicitRHS3d(:,:,7) + dt * IMa(intRK,1)*ImplicitRHS3d(:,:,4) + dt * IMa(intRK,2)*ImplicitRHS3d(:,:,5)+...
            dt * IMa(intRK,3)*ImplicitRHS3d(:,:,6);       
        %Information about the 2d mesh is contained in meshUnion
        [ImplicitRHS3d(:,:,intRK), ImplicitRHS3d(:,:,intRK+3), fphys{1}(:,:,1), fphys{1}(:,:,2)] = ...
            matUpdateImplicitVerticalDiffusion(SystemRHS, DiffusionCoefficient, obj.meshUnion(1), fphys{1});              
        
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
        
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
        % fphys2d = obj.matEvaluatePostFunc( fphys2d );
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
        %>Actually, boundary condition need to be imposed here
        %         obj.matUpdateExternalField( time + dt, fphys2d, fphys3d );    
        [ExplicitRHS2d(:,:,4), ExplicitRHS3d(:,:,4), ExplicitRHS3d(:,:,4+4)] = ...
            matCalculateExplicitRHSTerm(obj, fphys2d, fphys, obj.fext2d);    
        %>Update the velocity
        fphys{1}(:,:,1) = Tempfphys(:,:,1) + dt * EXb(1) * ExplicitRHS3d(:,:,1) + dt * EXb(2) * ExplicitRHS3d(:,:,2)+...
            dt * EXb(3) * ExplicitRHS3d(:,:,3) + dt * EXb(4) * ExplicitRHS3d(:,:,4) + dt * IMb(1) * ImplicitRHS3d(:,:,1)+...
            dt * IMb(2) * ImplicitRHS3d(:,:,2) + dt * IMb(3) * ImplicitRHS3d(:,:,3);
        
        fphys{1}(:,:,2) = Tempfphys(:,:,2) + dt * EXb(1) * ExplicitRHS3d(:,:,5) + dt * EXb(2) * ExplicitRHS3d(:,:,6)+...
            dt * EXb(3) * ExplicitRHS3d(:,:,7) + dt * EXb(4) * ExplicitRHS3d(:,:,8) + dt * IMb(1) * ImplicitRHS3d(:,:,4)+...
            dt * IMb(2) * ImplicitRHS3d(:,:,5) + dt * IMb(3) * ImplicitRHS3d(:,:,6);  
        
        fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXb(1) * ExplicitRHS2d(:,:,1) + dt * EXb(2) * ExplicitRHS2d(:,:,2)+...
            dt * EXb(3) * ExplicitRHS2d(:,:,3) + dt * EXb(4) * ExplicitRHS2d(:,:,4);        
        
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
    ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,4);
    ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 4*obj.Nvar);
    ImplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
    time = time + dt;
    
    fphys{1}(:,:,5) = obj.EddyViscositySolver.matUpdateEddyViscosity( obj, obj.mesh2d, obj.meshUnion(1), fphys2d, fphys, dt , time );
    %> Update the diffusion coefficient
    DiffusionCoefficient = fphys{1}(:,:,5)./fphys{1}(:,:,4).^2;
    display(time);
    obj.matUpdateOutputResult( time, fphys2d, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys2d = fphys2d;
obj.fphys = fphys;
obj.matUpdateFinalResult( time, fphys2d, fphys );
% obj.outputFile.closeOutputFile();
end

function [ExplicitRHS2d, ExplicitHuRHS3d, ExplicitHvRHS3d] = matCalculateExplicitRHSTerm(obj, fphys2d, fphys, fext2d)
  %> update the vertical velocity
    fphys{1}(:,:,3) = obj.matEvaluateVerticalVelocity( obj.meshUnion(1), fphys2d{1}, fphys{1} );
    obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d, fphys, fext2d);
    ExplicitRHS2d = obj.frhs2d{1}(:,:,1);
    obj.advectionSolver.evaluateAdvectionRHS( fphys );
    obj.matEvaluateSourceTerm( fphys );
    ExplicitHuRHS3d = obj.frhs{1}(:,:,1);
    ExplicitHvRHS3d = obj.frhs{1}(:,:,2);
end

function [ImplicithuRHS3d, ImplicithvRHS3d, hu, hv] = matUpdateImplicitVerticalDiffusion(SystemRHS, DiffusionCoefficient, meshUnion, fphys)


end

function [Explicita, Implicita, Explicitb, Implicitb, Parameterc] = GetRKParamter()
GAMA = 0.435866521508460;
beta1 = 1.369257448854393;
beta2 = -0.644363170684471;
alpha1 = -0.35;
alpha2 = -0.989175724679855;
Parameterc = [0 GAMA (1+GAMA)/2 1];
Explicita = [0 0 0 0;
    GAMA 0 0 0;
    (1+GAMA)/2-alpha1 alpha1 0 0;
    0 1-alpha2 alpha2 0];
Implicita = [GAMA 0 0;
    (1-GAMA)/2 GAMA 0;
    beta1 beta2 GAMA];
Explicitb = [0 beta1 beta2 GAMA];
Implicitb = [beta1 beta2 GAMA];
end