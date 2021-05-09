function matTimeSteppingNew( obj )
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys2d = obj.fphys2d;
fphys = obj.fphys;
visual = Visual2d( obj.mesh2d );
hwait = waitbar(0,'Runing MatSolver....');

% while( time < ftime )
%     dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
%     %       dt = 0.1;
%     if( time + dt > ftime )
%         dt = ftime - time;
%     end
%     
%     Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
%     Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
%     
%     tloc = time + dt;
%     obj.matUpdateExternalField( tloc, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     % $H^*u^*$, $H^*v^*$, $H^*w^*$
%     [ fphys{1}(:,:,obj.varFieldIndex)] = Tempfphys + dt * obj.frhs{1};
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
%     
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     % $H^*$
%     fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );    
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     [ fphys{1}(:,:,obj.varFieldIndex)] = fphys{1}(:,:,obj.varFieldIndex) + dt * obj.frhs{1};
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
%     %> The final velocity
%     [ fphys{1}(:,:,obj.varFieldIndex)] = 1/2 * Tempfphys + 1/2 * fphys{1}(:,:,obj.varFieldIndex);
%     
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     
%     fphys2d{1}(:,:,1) = fphys2d{1}(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
% %     fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
%     fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     %     visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
%     visual.drawResult( fphys2d{1}(:,:,1) );
% %     disp(max(max(fphys2d{1}(:,:,1))));
%     % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
%     %> reallocate the space for the rhs
%     time = time + dt;
%     
%     %> Update the diffusion coefficient
%     display(time);
%     obj.matUpdateOutputResult( time, fphys2d, fphys );
%     timeRatio = time / ftime;
%     waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
% end
% hwait.delete();
% obj.fphys2d = fphys2d;
% obj.fphys = fphys;
% obj.matUpdateFinalResult( time, fphys2d, fphys );
% obj.matClearGlobalMemory();
% % obj.outputFile.closeOutputFile();
% end



% while( time < ftime )
%     dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
%     %       dt = 0.1;
%     if( time + dt > ftime )
%         dt = ftime - time;
%     end
%     
%     Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
%     Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
%     
%     tloc = time + dt;
%     obj.matUpdateExternalField( tloc, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     % $H^*u^*$, $H^*v^*$, $H^*w^*$
%     [ fphys{1}(:,:,obj.varFieldIndex)] = Tempfphys + dt * obj.frhs{1};
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
%     
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
% 
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     % $H^*$
%     fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );    
%     
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     [ fphys{1}(:,:,obj.varFieldIndex)] = fphys{1}(:,:,obj.varFieldIndex) + dt * obj.frhs{1};
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
%     
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     
%     fphys2d{1}(:,:,1) = fphys2d{1}(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     %> The final velocity
%     [ fphys{1}(:,:,obj.varFieldIndex)] = 1/2 * Tempfphys + 1/2 * fphys{1}(:,:,obj.varFieldIndex);
%     
%     fphys2d{1}(:,:,1) = 1/2 * Tempfphys2d + 1/2 * fphys2d{1}(:,:,1);
%     
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
%     fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     %     visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
%     visual.drawResult( fphys2d{1}(:,:,1) );
%     disp(max(max(fphys2d{1}(:,:,1))));
%     % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
%     %> reallocate the space for the rhs
%     time = time + dt;
%     
%     %> Update the diffusion coefficient
%     %     display(time);
%     obj.matUpdateOutputResult( time, fphys2d, fphys );
%     timeRatio = time / ftime;
%     waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
% end
% hwait.delete();
% obj.fphys2d = fphys2d;
% obj.fphys = fphys;
% obj.matUpdateFinalResult( time, fphys2d, fphys );
% obj.matClearGlobalMemory();
% % obj.outputFile.closeOutputFile();
% end



% while( time < ftime )
%     dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
%     %       dt = 0.1;
%     if( time + dt > ftime )
%         dt = ftime - time;
%     end
%     
%     Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
%     Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
%     
%     tloc = time + dt;
%     obj.matUpdateExternalField( tloc, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     % $H^*u^*$, $H^*v^*$, $H^*w^*$
%     [ fphys{1}(:,:,obj.varFieldIndex)] = Tempfphys + dt * obj.frhs{1};
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     % $H^*$
%     fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );       
%     
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d ); 
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     [ fphys{1}(:,:,obj.varFieldIndex)] = fphys{1}(:,:,obj.varFieldIndex) + dt * obj.frhs{1};
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );    
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     
%     fphys2d{1}(:,:,1) = fphys2d{1}(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     %> The final velocity
%     [ fphys{1}(:,:,obj.varFieldIndex)] = 1/2 * Tempfphys + 1/2 * fphys{1}(:,:,obj.varFieldIndex);
%     
%     fphys2d{1}(:,:,1) = 1/2 * Tempfphys2d + 1/2 * fphys2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
%     fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);    
%     
%     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     %     visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
%     visual.drawResult( fphys2d{1}(:,:,1) );
% %     disp(max(max(fphys2d{1}(:,:,1))));
%     % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
%     %> reallocate the space for the rhs
%     time = time + dt;
%     
%     %> Update the diffusion coefficient
%     display(time);
%     obj.matUpdateOutputResult( time, fphys2d, fphys );
%     timeRatio = time / ftime;
%     waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
% end
% hwait.delete();
% obj.fphys2d = fphys2d;
% obj.fphys = fphys;
% obj.matUpdateFinalResult( time, fphys2d, fphys );
% obj.matClearGlobalMemory();
% % obj.outputFile.closeOutputFile();
% end



%%下面的时间递进格式可以算的稳(100S)，但是耗散大的无可理喻
% while( time < ftime )
%     dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
%     %       dt = 0.1;
%     if( time + dt > ftime )
%         dt = ftime - time;
%     end
%     
%     Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
%     Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
%     
%     tloc = time + dt;
%     obj.matUpdateExternalField( tloc, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     % $H^*u^*$, $H^*v^*$, $H^*w^*$
%     [ fphys{1}(:,:,obj.varFieldIndex)] = Tempfphys + dt * obj.frhs{1};
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     % $H^*$
%     fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) ); 
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
%     
% %     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d ); 
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     [ fphys{1}(:,:,obj.varFieldIndex)] = fphys{1}(:,:,obj.varFieldIndex) + dt * obj.frhs{1}; 
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     
%     fphys2d{1}(:,:,1) = fphys2d{1}(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
%     
%     fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );   
%     %> The final velocity
%     [ fphys{1}(:,:,obj.varFieldIndex)] = 1/2 * Tempfphys + 1/2 * fphys{1}(:,:,obj.varFieldIndex);
%     
%     fphys2d{1}(:,:,1) = 1/2 * Tempfphys2d + 1/2 * fphys2d{1}(:,:,1); 
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
%     
%     fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
%     
% %     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     %     visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
%     visual.drawResult( fphys2d{1}(:,:,1) );
% %     disp(max(max(fphys2d{1}(:,:,1))));
%     % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
%     %> reallocate the space for the rhs
%     time = time + dt;
%     
%     %> Update the diffusion coefficient
%     display(time);
%     obj.matUpdateOutputResult( time, fphys2d, fphys );
%     timeRatio = time / ftime;
%     waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
% end
% hwait.delete();
% obj.fphys2d = fphys2d;
% obj.fphys = fphys;
% obj.matUpdateFinalResult( time, fphys2d, fphys );
% obj.matClearGlobalMemory();
% % obj.outputFile.closeOutputFile();
% end




%%20s测试算的稳，但是耗散还是很大很大
% while( time < ftime )
%     dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
%     %       dt = 0.1;
%     if( time + dt > ftime )
%         dt = ftime - time;
%     end
%     
%     Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
%     Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
%     
%     tloc = time + dt;
%     obj.matUpdateExternalField( tloc, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     % $H^*u^*$, $H^*v^*$, $H^*w^*$
%     [ fphys{1}(:,:,obj.varFieldIndex)] = Tempfphys + dt * obj.frhs{1};
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     % $H^*$
%     fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) ); 
%     
%     fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );    
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%     
% %     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d ); 
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
%     
%     obj.matEvaluateSourceTerm( fphys );
%     
%     [ fphys{1}(:,:,obj.varFieldIndex)] = fphys{1}(:,:,obj.varFieldIndex) + dt * obj.frhs{1}; 
%     
%     obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
%     
%     fphys2d{1}(:,:,1) = fphys2d{1}(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
%     
%     fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6); 
%     
%     fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt ); 
%     
%     fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
%     fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
%       
%     %> The final velocity
%     [ fphys{1}(:,:,obj.varFieldIndex)] = 1/2 * Tempfphys + 1/2 * fphys{1}(:,:,obj.varFieldIndex);
%     
%     fphys2d{1}(:,:,1) = 1/2 * Tempfphys2d + 1/2 * fphys2d{1}(:,:,1); 
%     
%     fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
%     
%     fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
%     
% %     fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
%     
%     %> update the vertical velocity
%     fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
%     
%     %     visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
%     visual.drawResult( fphys2d{1}(:,:,1) );
% %     disp(max(max(fphys2d{1}(:,:,1))));
%     % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
%     %> reallocate the space for the rhs
%     time = time + dt;
%     
%     %> Update the diffusion coefficient
%     display(time);
%     obj.matUpdateOutputResult( time, fphys2d, fphys );
%     timeRatio = time / ftime;
%     waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
% end
% hwait.delete();
% obj.fphys2d = fphys2d;
% obj.fphys = fphys;
% obj.matUpdateFinalResult( time, fphys2d, fphys );
% obj.matClearGlobalMemory();
% % obj.outputFile.closeOutputFile();
% end


%% 20s 测试模型稳定，一层模型垂向三阶近似，计算结果很好，耗散较汇报更小。
while( time < ftime )
    dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
    %       dt = 0.1;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
    Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
    
    tloc = time + dt;
    obj.matUpdateExternalField( tloc, fphys2d, fphys );
    
    obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
    
    obj.matEvaluateSourceTerm( fphys );
    
    % $H^*u^*$, $H^*v^*$, $H^*w^*$
    [ fphys{1}(:,:,obj.varFieldIndex)] = Tempfphys + dt * obj.frhs{1};
    
    obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
    % $H^*$
    fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) ); 
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
    
    fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );    
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    
    fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d ); 
    
    %> update the vertical velocity
    fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    
    obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
    
    obj.matEvaluateSourceTerm( fphys );
    
    [ fphys{1}(:,:,obj.varFieldIndex)] = fphys{1}(:,:,obj.varFieldIndex) + dt * obj.frhs{1}; 
    
    obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
    
    fphys2d{1}(:,:,1) = fphys2d{1}(:,:,1) + dt * obj.frhs2d{1}(:,:,1);
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6); 
    
    fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata( obj, fphys, fphys2d, dt );
    
    fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
      
    %> The final velocity
    [ fphys{1}(:,:,obj.varFieldIndex)] = 1/2 * Tempfphys + 1/2 * fphys{1}(:,:,obj.varFieldIndex);
    
    fphys2d{1}(:,:,1) = 1/2 * Tempfphys2d + 1/2 * fphys2d{1}(:,:,1); 
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);  
    
    fphys = obj.NonhydrostaticSolver.matUpdataVerticalVelocity( obj, fphys, fphys2d );
    
    %> update the vertical velocity
    fphys{1}(:,:,3) = obj.VerticalVelocitySolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
    
    %     visual.drawResult( fphys2d{1}(:,:,1) + fphys2d{1}(:,:,4) );
    visual.drawResult( fphys2d{1}(:,:,1) );
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