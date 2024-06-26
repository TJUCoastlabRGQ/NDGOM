function matEvaluateSSPRK22(obj)
profile on;
[ rkb, rkt] = GetRKParameter();
fphys2d = obj.fphys2d;
fphys = obj.fphys;
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
obj.ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,2);
obj.ExplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 2 * obj.Nvar);
%Actually, the following variable is not used, it is allocated here to
% accomadate the usage of other time stepping method such as IMEXRK222, note added on 20211223
obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nvar);
visual = Visual2d( obj.mesh2d );
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = 0.3 * obj.matUpdateTimeInterval( fphys2d );
%     dt = 1.5;
%     dt = 10;
    %       dt = 0.1;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
    Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
    
    for intRK = 1:1
        tloc = time + rkt(intRK) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys2d, fphys );
        
        obj.matCalculateExplicitRHSTerm( fphys2d, fphys, 2, intRK);
        
%         fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + rkb(intRK, 1) * dt * obj.ExplicitRHS2d(:,:,1) + rkb(intRK, 2) * dt * obj.ExplicitRHS2d(:,:,2);
%         
%         fphys{1}(:,:,obj.varFieldIndex) = Tempfphys + rkb(intRK, 1) * dt * obj.ExplicitRHS(:,:,1:obj.Nvar:end) + rkb(intRK, 2) * dt * obj.ExplicitRHS(:,:,2:obj.Nvar:end);

        fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) +  dt * obj.ExplicitRHS2d(:,:,1);
        
        fphys{1}(:,:,obj.varFieldIndex) = Tempfphys +  dt * obj.ExplicitRHS(:,:,1:2:end);
        
%         [ fphys ] = obj.matImposeLimiter( fphys );
        
%         [ fphys ] = obj.limiter.matLimitNew(obj, fphys);
        
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
        
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
        
        fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
        %         [ fphys ] = obj.matFilterSolution( fphys );
        
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    
    for intRK = 2:2
        tloc = time + rkt(intRK) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys2d, fphys );
        
        obj.matCalculateExplicitRHSTerm( fphys2d, fphys, 2, intRK);
        
        fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + rkb(intRK, 1) * dt * obj.ExplicitRHS2d(:,:,1) + rkb(intRK, 2) * dt * obj.ExplicitRHS2d(:,:,2);
        
        fphys{1}(:,:,obj.varFieldIndex) = Tempfphys + rkb(intRK, 1) * dt * obj.ExplicitRHS(:,:,1:2:end) + rkb(intRK, 2) * dt * obj.ExplicitRHS(:,:,2:2:end);
        
%         [ fphys ] = obj.matImposeLimiter( fphys );
%         [ fphys ] = obj.limiter.matLimitNew(obj, fphys);
        
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
        %         [ fphys ] = obj.matFilterSolution( fphys );
    end
    
    [ fphys{1}(:,:,obj.varFieldIndex)] = ...
        obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
        fphys2d{1}(:,:,1), fphys{1}(:,:,4), fphys{1}(:,:,obj.varFieldIndex), 1, dt, 1,...
        2, fphys{1}(:,:,1), fphys{1}(:,:,2), time, fphys );
    
%     [ fphys{1}(:,:,5)] = obj.VerticalEddyViscositySolver.nv;
%     
%     [ fphys{1}(:,:,14)] = obj.VerticalEddyViscositySolver.Tke;
%     
%     [ fphys{1}(:,:,15)] = obj.VerticalEddyViscositySolver.Eps;
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
    fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    
    visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys2d, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys2d = fphys2d;
obj.fphys = fphys;
obj.matUpdateFinalResult( time, fphys2d, fphys );
profile viewer;
end

function [ rkb, rkt] = GetRKParameter
rkb = [1 0;
    1/2 1/2];
rkt = [0 1];
end