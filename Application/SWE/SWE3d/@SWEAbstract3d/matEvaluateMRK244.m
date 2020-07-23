function matEvaluateMRK244( obj )
c = [1/4 1/3 1/2 1];
Stage = size(c,2);
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys2d = obj.fphys2d;
fphys = obj.fphys;
%> allocate space for the rhs to be stored
obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, Stage  * obj.Nvar);
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
    
    for intRK = 1:4
        tloc = time + c( intRK ) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys2d, fphys );
        
        obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );

        fphys2d{1}(:,:,1) = Tempfphys2d + c(intRK) * dt * obj.frhs2d{1};
        
        obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
        
        SystemRHS = Tempfphys + c( intRK ) * dt * obj.frhs{1};
        
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
        
        %dt here is problematic
        [ fphys{1}(:,:,obj.varFieldIndex)] = ...
            obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
            fphys2d{1}(:,:,1), fphys{1}(:,:,4), SystemRHS, c(intRK), dt, intRK,...
            Stage, fphys{1}(:,:,1), fphys{1}(:,:,2), time );
        
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
        
        
        %> update the vertical velocity
        fphys{1}(:,:,3) = obj.VertSolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );

                
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
        %         fphys2d = obj.matEvaluatePostFunc( fphys2d );
        
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    visual.drawResult( fphys2d{1}(:,:,1) );
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    time = time + dt;
    
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