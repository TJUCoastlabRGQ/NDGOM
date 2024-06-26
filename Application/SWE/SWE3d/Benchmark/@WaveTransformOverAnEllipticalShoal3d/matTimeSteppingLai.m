function matTimeSteppingLai( obj )
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys2d = obj.fphys2d;
fphys = obj.fphys;
visual = Visual2d( obj.mesh2d );
hwait = waitbar(0,'Runing MatSolver....');
[time, fphys2d{1}, fphys{1}] = obj.NonhydrostaticSolver.matHotStart(obj.meshUnion(1),'C:\Students\RGQ\FinalCaseResult\2d\WaveTransformOverAnEllipticalShoal3d.1-1.1000.nc',...
'C:\Students\RGQ\FinalCaseResult\3d\WaveTransformOverAnEllipticalShoal3d.1-1.1000.nc', fphys2d{1}, fphys{1});
while( time < ftime )
    dt = 0.5*obj.matUpdateTimeInterval( fphys );
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
    
    obj.matEvaluateSourceTerm( fphys );
    
    obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d, fphys );
    
    % $H^*u^*$, $H^*v^*$, $H^*w^*$
    [ fphys{1}(:,:,obj.varFieldIndex)] = Tempfphys + dt * obj.frhs{1};
    
    % $H^*$
    fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) ); 
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
    
    fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt ); 
    
%     [ fphys ] = obj.Limiter.matLimitNew( obj, fphys );  
    
    %进去以后计算底部流速
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d ); 
    
    %> update the vertical velocity
    fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    
    obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
    
    obj.matEvaluateSourceTerm( fphys );
    
    obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d, fphys );
    
    [ fphys{1}(:,:,obj.varFieldIndex)] = fphys{1}(:,:,obj.varFieldIndex) + dt * obj.frhs{1}; 
    
    fphys2d{1}(:,:,1) = fphys2d{1}(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6); 
    
    fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
    
%     [ fphys ] = obj.Limiter.matLimitNew( obj, fphys );
    
%     [ fphys ] = obj.matImposeLimiter( fphys );  
    
    %进去以后计算底部流速，并利用上一个中间步的流速进行底部边界施加边界条件
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
      
    %> The final velocity
    [ fphys{1}(:,:,obj.varFieldIndex)] = 1/2 * Tempfphys + 1/2 * fphys{1}(:,:,obj.varFieldIndex);
    
%     [ fphys ] = obj.Limiter.matLimitNew( obj, fphys );
    
    fphys2d{1}(:,:,1) = 1/2 * Tempfphys2d + 1/2 * fphys2d{1}(:,:,1); 
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
    
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
    
    %> update the vertical velocity
    fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    
    visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
%     visual.drawResult( fphys2d{1}(:,:,1) );
%     disp(max(max(fphys2d{1}(:,:,1))));
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
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
obj.matClearGlobalMemory();
% obj.outputFile.closeOutputFile();
end