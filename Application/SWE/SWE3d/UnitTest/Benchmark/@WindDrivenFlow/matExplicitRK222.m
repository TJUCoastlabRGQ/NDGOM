function matExplicitRK222( obj )
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
time = obj.startTime;
ftime = obj.finalTime;

fphys2d = obj.fphys2d;
fphys = obj.fphys;
%> allocate space for the rhs to be stored
ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,3);
ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
visual = Visual2d( obj.mesh2d );
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    dt = 1/100 * 0.4 * obj.matUpdateTimeInterval( fphys2d );
%       dt = 0.1;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,1);
    Tempfphys = fphys{1}(:,:,1:2);
    
    tloc = time + c( 1 ) * dt;
    obj.matUpdateExternalField( tloc, fphys2d, fphys );
    [ExplicitRHS2d(:,:,1), ExplicitRHS3d(:,:,1), ExplicitRHS3d(:,:,1+3)] = ...
        matCalculateExplicitRHSTerm(obj, fphys2d, fphys, obj.fext2d);
    for intRK = 1:2
        tloc = time + c( intRK+1 ) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys2d, fphys );
        
        fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXa(intRK+1,1) * ExplicitRHS2d(:,:,1) + dt * EXa(intRK+1,2) * ExplicitRHS2d(:,:,2);
        fphys{1}(:, :, 1) = Tempfphys(:,:,1) + dt * EXa(intRK+1,1) * ExplicitRHS3d(:,:,1) + dt * EXa(intRK+1,2) * ExplicitRHS3d(:,:,2);
        fphys{1}(:, :, 2) = Tempfphys(:,:,2) + dt * EXa(intRK+1,1) * ExplicitRHS3d(:,:,4) + dt * EXa(intRK+1,2) * ExplicitRHS3d(:,:,5);
        
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
        
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
        
        %> update the vertical velocity
        fphys{1}(:,:,3) = obj.matEvaluateVerticalVelocity( obj.meshUnion(1), fphys2d, fphys );
        
        [ExplicitRHS2d(:,:,intRK+1), ExplicitRHS3d(:,:,intRK+1), ExplicitRHS3d(:,:,intRK+1+3)] = ...
            matCalculateExplicitRHSTerm(obj, fphys2d, fphys, obj.fext2d);
        
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
%         fphys2d = obj.matEvaluatePostFunc( fphys2d );
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    %>Update the velocity
    fphys{1}(:,:,1) = Tempfphys(:,:,1) + dt * EXb(1) * ExplicitRHS3d(:,:,1) + dt * EXb(2) * ExplicitRHS3d(:,:,2)+...
        dt * EXb(3) * ExplicitRHS3d(:,:,3);
    
    fphys{1}(:,:,2) = Tempfphys(:,:,2) + dt * EXb(1) * ExplicitRHS3d(:,:,4) + dt * EXb(2) * ExplicitRHS3d(:,:,5)+...
        dt * EXb(3) * ExplicitRHS3d(:,:,6);
    
    fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXb(1) * ExplicitRHS2d(:,:,1) + dt * EXb(2) * ExplicitRHS2d(:,:,2)+...
        dt * EXb(3) * ExplicitRHS2d(:,:,3);
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
    fphys{1}(:,:,3) = obj.matEvaluateVerticalVelocity( obj.meshUnion(1), fphys2d, fphys );
    visual.drawResult( fphys2d{1}(:,:,1) );
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
    ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,3);
    ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
    time = time + dt;
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
    %     [fphys{1}(:,:,5), obj.Cf{1}] = obj.EddyViscositySolver.matUpdateEddyViscosity( obj, obj.mesh2d, ...
    %         obj.meshUnion(1), fphys2d, fphys, dt , time, obj.WindTaux{1}, obj.WindTauy{1} );
    
    %> Update the diffusion coefficient
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
obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
obj.viscositySolver.matEvaluateRHS( fphys );
obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
obj.matEvaluateSourceTerm( fphys );
ExplicitRHS2d = obj.frhs2d{1}(:,:,1);
ExplicitHuRHS3d = obj.frhs{1}(:,:,1);
ExplicitHvRHS3d = obj.frhs{1}(:,:,2);
end


function [Explicita, Implicita, Explicitb, Implicitb, Parameterc] = GetRKParamter()
GAMA = (2-sqrt(2))/2;
delta = 1-1/(2*GAMA);
Parameterc = [0 GAMA 1];
Explicita = [0 0 0;
    GAMA 0 0;
    delta 1-delta 0];
Implicita = [GAMA 0;
    (1-GAMA) GAMA];
Explicitb = [delta 1-delta 0];
Implicitb = [1-GAMA GAMA];
end