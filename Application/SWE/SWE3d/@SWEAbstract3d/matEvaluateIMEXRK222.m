function matEvaluateIMEXRK222( obj )
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
Stage = size(EXa,2);
time = obj.startTime;
ftime = obj.finalTime;
fphys2d = obj.fphys2d;
fphys = obj.fphys;
%> allocate space for the rhs to be stored
obj.ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,Stage);
obj.ExplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, Stage*obj.Nvar);
obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, ( Stage - 1 ) * obj.Nvar);
SystemRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nvar);
visual = Visual2d( obj.mesh2d );
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
    %       dt = 0.1;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
    Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
    
    tloc = time + c( 1 ) * dt;
    obj.matUpdateExternalField( tloc, fphys2d, fphys );
    obj.matCalculateExplicitRHSTerm( fphys2d, fphys, Stage, 1);
    for intRK = 1:2
        tloc = time + c( intRK+1 ) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys2d, fphys );
        
        SystemRHS = obj.matAssembleSystemRHS( Tempfphys, SystemRHS, EXa(intRK+1,:), IMa(intRK,:), dt);
        
        fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXa(intRK+1,1) * obj.ExplicitRHS2d(:,:,1) +...
            dt * EXa(intRK+1,2) * obj.ExplicitRHS2d(:,:,2) + dt * EXa(intRK+1,3) * obj.ExplicitRHS2d(:,:,3);   
        
        %dt here is problematic
        [ fphys{1}(:,:,obj.varFieldIndex)] = ...
            obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
            fphys2d{1}(:,:,1), fphys{1}(:,:,4), SystemRHS, IMa(intRK,intRK), dt, intRK,...
            Stage, fphys{1}(:,:,1), fphys{1}(:,:,2), time );
        
%         [ fphys ] = obj.matImposeLimiter( fphys );  
%         [ fphys ] = obj.matFilterSolution( fphys ); 
        
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
        
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
        
        %> update the vertical velocity
        fphys{1}(:,:,3) = obj.VertSolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );

        
        obj.matCalculateExplicitRHSTerm( fphys2d, fphys, Stage, intRK + 1);
        
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
        %         fphys2d = obj.matEvaluatePostFunc( fphys2d );
        
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    %>Update the velocity
    
    for i = 1:obj.Nvar
        fphys{1}(:,:,obj.varFieldIndex(i)) = Tempfphys(:,:,i) + dt * EXb(1) * obj.ExplicitRHS(:,:,(i-1)*3+1) + dt * EXb(2) * obj.ExplicitRHS(:,:,(i-1)*3+2)+...
            dt * EXb(3) * obj.ExplicitRHS(:,:,(i-1)*3+3) + dt * IMb(1) * obj.ImplicitRHS(:,:,(i-1)*2+1) + dt * IMb(2) * obj.ImplicitRHS(:,:,(i-1)*2+2);
        %         if(any(abs(obj.ImplicitRHS3d(:,:,(i-1)*2+1)) > 1e-10))
        %             fprintf('Nonzero element contained\n');
        %         elseif(any(abs(obj.ImplicitRHS3d(:,:,(i-1)*2+2)) > 1e-10))
        %             fprintf('Nonzero element contained\n');
        %         end
    end
    
    fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXb(1) * obj.ExplicitRHS2d(:,:,1) + dt * EXb(2) * obj.ExplicitRHS2d(:,:,2)+...
        dt * EXb(3) * obj.ExplicitRHS2d(:,:,3);
    
%    [ fphys ] = obj.matImposeLimiter( fphys );   
%     [ fphys ] = obj.matFilterSolution( fphys ); 
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
    visual.drawResult( fphys2d{1}(:,:,1) );
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
    obj.ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,Stage);
    obj.ExplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, Stage*obj.Nvar);
    obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, ( Stage - 1 ) * obj.Nvar);
    time = time + dt;
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
    fphys{1}(:,:,3) = obj.VertSolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    
    %> Update the diffusion coefficient
%     display(time);
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

function fphys = matFilterSolution( obj, fphys )
Coef = diag([1,0.9,0.8,0.7,0.6,0.5,0.4,0.3]);
FilterMatrix = (obj.meshUnion.cell.V * Coef) / obj.meshUnion.cell.V;
fphys{1}(:,:,1) = FilterMatrix * fphys{1}(:,:,1);
fphys{1}(:,:,2) = FilterMatrix * fphys{1}(:,:,2);
end

