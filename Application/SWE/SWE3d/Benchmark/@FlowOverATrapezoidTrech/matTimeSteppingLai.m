function matTimeSteppingLai( obj )
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys2d = obj.fphys2d;
fphys = obj.fphys;
visual = Visual2d( obj.mesh2d );
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = 0.3 * obj.matUpdateTimeInterval( fphys );
%     dt = 0.005;
    %       dt = 0.1;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
    Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
    
   obj.NonhydrostaticSolver.matCalculateBottomVerticalVelocity( obj, fphys );

    
    tloc = time + dt;
    obj.matUpdateExternalField( tloc, fphys2d, fphys );
    
    obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
    
    obj.HorizontalEddyViscositySolver.matEvaluateDiffRHS( obj, fphys{1} );
    
    obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d, fphys );
    
    obj.matEvaluateSourceTerm( fphys );
    
    % $H^*u^*$, $H^*v^*$, $H^*w^*$
    [ fphys{1}(:,:,obj.varFieldIndex)] = Tempfphys + dt * obj.frhs{1}; 
    
    disp("Fot the right hand side of hu\n")
    DataCheck(obj.frhs{1}(:,:,1), obj.meshUnion.K/2);
    
    disp("Fot the right hand side of hv\n")
    DataCheck(obj.frhs{1}(:,:,2), obj.meshUnion.K/2);
    
    % $H^*$
    fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
    
    disp("Fot the right hand side of h\n")
    DataCheck(obj.frhs2d{1}(:,:,1), obj.meshUnion.mesh2d.K/2);
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) ); 
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
    
    DataCheck(fphys{1}(: , :, 7), obj.meshUnion.K/2);
    
    [ fphys ] = obj.matImposeLimiter( fphys );
    
    [ fphys{1}(:,:,obj.varFieldIndex)] = ...
        obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
        fphys2d{1}(:,:,1), fphys{1}(:,:,4), fphys{1}(:,:,obj.varFieldIndex), 1, dt, 1,...
        2, fphys{1}(:,:,1), fphys{1}(:,:,2), time, fphys );    
    
%     obj.NonhydrostaticSolver.matCalculateNonhydroRHS(obj, fphys, fphys2d );
    
     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );    
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d ); 
    
    %> update the vertical velocity
    fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    
    obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
    
    obj.HorizontalEddyViscositySolver.matEvaluateDiffRHS( obj, fphys{1} );
    
    obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d, fphys );
    
    obj.matEvaluateSourceTerm( fphys );
    
    [ fphys{1}(:,:,obj.varFieldIndex)] = fphys{1}(:,:,obj.varFieldIndex) + dt * obj.frhs{1}; 
    
    disp("Fot the right hand side of hu\n")
    DataCheck(obj.frhs{1}(:,:,1), obj.meshUnion.K/2);
    disp("Fot the right hand side of hv\n")
    DataCheck(obj.frhs{1}(:,:,2), obj.meshUnion.K/2);
    
    fphys2d{1}(:,:,1) = fphys2d{1}(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
    disp("Fot the right hand side of h\n")
    DataCheck(obj.frhs2d{1}(:,:,1), obj.meshUnion.mesh2d.K/2);
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
    
    DataCheck(fphys{1}(: , :, 7), obj.meshUnion.K/2);
    
    [ fphys ] = obj.matImposeLimiter( fphys );
    
    [ fphys{1}(:,:,obj.varFieldIndex)] = ...
        obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
        fphys2d{1}(:,:,1), fphys{1}(:,:,4), fphys{1}(:,:,obj.varFieldIndex), 1, dt, 1,...
        2, fphys{1}(:,:,1), fphys{1}(:,:,2), time, fphys );       
    
%     obj.NonhydrostaticSolver.matCalculateNonhydroRHS(obj, fphys, fphys2d );
    
    fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
    
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
      
    %> The final velocity
    [ fphys{1}(:,:,obj.varFieldIndex)] = 1/2 * Tempfphys + 1/2 * fphys{1}(:,:,obj.varFieldIndex);
    
    fphys2d{1}(:,:,1) = 1/2 * Tempfphys2d + 1/2 * fphys2d{1}(:,:,1); 
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
    
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
    
    %> update the vertical velocity
    fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    
%     visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
    visual.drawResult( fphys2d{1}(:,:,1) );
%     disp(max(max(fphys2d{1}(:,:,1))));
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
    time = time + dt;
    
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
obj.matClearGlobalMemory();
% obj.outputFile.closeOutputFile();
end