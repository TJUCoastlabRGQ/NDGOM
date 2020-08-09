function matEvaluateIMEXRK222( obj )
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
Stage = size(EXa,2);
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
dt = obj.getOption('timeInterval');
fphys = obj.fphys;
%> allocate space for the rhs to be stored
obj.ExplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, Stage*obj.Nvar);
obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, ( Stage - 1 ) * obj.Nvar);
SystemRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nvar);
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    %     dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
    %       dt = 0.1;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
    
    tloc = time + c( 1 ) * dt;
    obj.matUpdateExternalField( tloc, fphys );
    obj.matCalculateExplicitRHSTerm( fphys, tloc,  Stage, 1);
    for intRK = 1:2
        tloc = time + c( intRK+1 ) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys );
        
        SystemRHS = obj.matAssembleSystemRHS( Tempfphys, SystemRHS, EXa(intRK+1,:), IMa(intRK,:), dt);
        %da gai
        [ fphys{1}(:,:,obj.varFieldIndex)] = ...
            obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
            SystemRHS, IMa(intRK,intRK), dt, intRK, Stage );    
        
        obj.matCalculateExplicitRHSTerm( fphys, tloc , Stage, intRK + 1);
        
    end
    %>Update the velocity
    
    for i = 1:obj.Nvar
        fphys{1}(:,:,obj.varFieldIndex(i)) = Tempfphys(:,:,i) + dt * EXb(1) * obj.ExplicitRHS(:,:,(i-1)*3+1) + dt * EXb(2) * obj.ExplicitRHS(:,:,(i-1)*3+2)+...
            dt * EXb(3) * obj.ExplicitRHS(:,:,(i-1)*3+3) + dt * IMb(1) * obj.ImplicitRHS(:,:,(i-1)*2+1) + dt * IMb(2) * obj.ImplicitRHS(:,:,(i-1)*2+2);
    end
    
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
    obj.ExplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, Stage*obj.Nvar);
    obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, ( Stage - 1 ) * obj.Nvar);
    time = time + dt;
%     fphys{1}(:,:,1) = obj.matGetExtFunc( time );    
    %> Update the diffusion coefficient
    display(time);
    obj.matUpdateOutputResult( time,  fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys = fphys;
obj.matUpdateFinalResult( time, fphys );
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