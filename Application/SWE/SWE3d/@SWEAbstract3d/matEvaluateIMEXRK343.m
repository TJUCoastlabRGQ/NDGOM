function matEvaluateIMEXRK343( obj )
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');

fphys2d = obj.fphys2d;
fphys = obj.fphys;
%> allocate space for the rhs to be stored
obj.ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,4);
obj.ExplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, Stage*obj.Nvar);
obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, ( Stage - 1 ) * obj.Nvar);
SystemRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nvar);
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
    %     dt = 2.7;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,1);
    Tempfphys = fphys{1}(:,:,1:2);
    
    tloc = time + c( 1 ) * dt;
    obj.matUpdateExternalField( tloc, fphys2d, fphys );
    obj.matCalculateExplicitRHSTerm( fphys2d,  fphys, tloc,  Stage, 1);
    
    for intRK = 1:3
        tloc = time + c( intRK + 1 ) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys2d, fphys );
        %> Calculate the intermediate water depth
        SystemRHS = obj.matAssembleSystemRHS( Tempfphys, SystemRHS, EXa(intRK+1,:), IMa(intRK,:), dt);
        
        
        %> Calculate the right hand side for the global system about the three-dimensional horizontal momentum
        [ fphys{1}(:,:,obj.varFieldIndex)] = ...
            obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
            SystemRHS, IMa(intRK,intRK), dt, intRK, Stage );
        
        
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
        %> update the vertical velocity
        
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
        fphys{1}(:,:,3) = obj.VertSolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
        %> Calculation of the right hand side corresponds to the discretization of the non-stiff term at stage intRK+1 with the
        %> newly calculated intermediate water depth and three-dimensional horizontal momentum
        obj.matCalculateExplicitRHSTerm( fphys2d,  fphys, tloc,  Stage, intRK + 1);
        
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
        fphys2d = obj.matEvaluatePostFunc( fphys2d );
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    %> Update the water depth and the three-dimensional horizontal momentum at the next step
    for i = 1:obj.Nvar
        fphys{1}(:,:,obj.varFieldIndex(i)) = Tempfphys(:,:,i) + dt * EXb(1) * obj.ExplicitRHS(:,:,(i-1)*4+1) + dt * EXb(2) * obj.ExplicitRHS(:,:,(i-1)*4+2)+...
            dt * EXb(3) * obj.ExplicitRHS(:,:,(i-1)*4+3) + dt * EXb(4) * obj.ExplicitRHS(:,:,(i-1)*4+4) + ...
            dt * IMb(1) * obj.ImplicitRHS(:,:,(i-1)*3+1) + dt * IMb(2) * obj.ImplicitRHS(:,:,(i-1)*3+2) + dt * IMb(3) * obj.ImplicitRHS(:,:,(i-1)*3+3);
    end
    
    fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXb(1) * obj.ExplicitRHS2d(:,:,1) + dt * EXb(2) * obj.ExplicitRHS2d(:,:,2)+...
        dt * EXb(3) * obj.ExplicitRHS2d(:,:,3) + dt * EXb(4) * obj.ExplicitRHS2d(:,:,4);
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
    %> update the vertical velocity
    fphys{1}(:,:,3) = obj.VertSolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    %>Update the velocity
    visual.drawResult( fphys2d{1}(:,:,1) );
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
    obj.ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,4);
    obj.ExplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, Stage*obj.Nvar);
    obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, ( Stage - 1 ) * obj.Nvar);
    time = time + dt;
    
    %> Update the diffusion coefficient
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

function [Explicita, Implicita, Explicitb, Implicitb, Parameterc] = GetRKParamter()
data = roots([6 -18 9 -1]);
GAMA = data(2);
beta1 = -1.5*GAMA^2+4*GAMA-1/4;
beta2 = 1.5*GAMA^2-5*GAMA+5/4;
alpha1 = -0.35;
alpha2 = (1/3-2*GAMA^2-2*beta2*alpha1*GAMA)/(GAMA*(1-GAMA));
% GAMA = 0.435866521508460;
% beta1 = 1.208496649176012;
% beta2 = -0.644363170684471;
% alpha1 = -0.35;
% alpha2 = -0.989175724679855;
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