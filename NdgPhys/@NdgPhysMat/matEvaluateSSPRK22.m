function matEvaluateSSPRK22( obj )

tic; 
Nmesh = obj.Nmesh;

time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');

fphys0 = cell( obj.Nmesh, 1 );
for n = 1:obj.Nmesh
    fphys0{n} = zeros( obj.meshUnion(n).cell.Np, obj.meshUnion(n).K, obj.Nvar );
end

fphys = obj.fphys;
% init limiter and output file
visual = makeVisualizationFromNdgPhys( obj );
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = obj.matUpdateTimeInterval( fphys ) * 0.3;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for n = 1:obj.Nmesh
        fphys0{n} = fphys{n};
    end
    
    for intRK = 1:2
        tloc = time + dt;
        obj.matUpdateExternalField( tloc, fphys );
        obj.matEvaluateRHS( fphys );
        
        for n = 1:Nmesh
            %             resQ{n} = rk4a(intRK)*resQ{n} + dt*obj.frhs{n};
            fphys{n}(:,:, obj.varFieldIndex) ...
                = fphys{n}(:,:, obj.varFieldIndex) + dt*obj.frhs{n};
        end
        
%         fphys = obj.matEvaluateLimiter( fphys );
        fphys = obj.matEvaluatePostFunc( fphys );
        fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata(obj, fphys, dt);
    end
    
    for n = 1:obj.Nmesh
        fphys{n}(:,:, obj.varFieldIndex) = 1/2 * fphys0{n}(:,:, obj.varFieldIndex)...
            + 1/2 * fphys{n}(:,:, obj.varFieldIndex) ;
    end
    %     obj.meshUnion(1).draw( fphys{1}(:,:,1) );
    %     drawnow;
%     visual.drawResult( fphys{1}(:, :, 1) + fphys{1}(:, :, 4) );
        visual.drawResult( fphys{1}(:, :, 1) );
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ...
        ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
% hwait.delete();
obj.matUpdateFinalResult( time, fphys );
obj.fphys = fphys;
toc;
end

