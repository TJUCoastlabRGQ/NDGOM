function matEvaluateEuler( obj )

Nmesh = obj.Nmesh;

time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');

fphys = obj.fphys;

visual = makeVisualizationFromNdgPhys( obj );
% init limiter and output file
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = obj.matUpdateTimeInterval( fphys );
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    tloc = time + dt;
    obj.matUpdateExternalField( tloc, fphys );
    obj.matEvaluateRHS( fphys );
    
    for n = 1:Nmesh
        fphys{n}(:,:, obj.varFieldIndex) ...
            = fphys{n}(:,:, obj.varFieldIndex) + dt*obj.frhs{n};
    end
    
%     fphys = obj.matEvaluateLimiter( fphys );
    fphys = obj.matEvaluatePostFunc( fphys );
    
%     for m = 1:obj.Nmesh
%         obj.meshUnion(m).draw( fphys{m}(:,:,1) );
%     end
%     drawnow;
    
    time = time + dt;
%     obj.dt = dt;
    fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata(obj, fphys, dt);
    obj.matUpdateOutputResult( time, fphys );
    visual.drawResult( fphys{1}(:, :, 1) );
    waitbar( time/ftime, hwait, 'Runing MatSolver....');
%     obj.meshUnion(1).draw( fphys{1}(:,:,1) );
%     drawnow;    
end
hwait.delete();
obj.matUpdateFinalResult( time, fphys );
obj.fphys = fphys;
end

