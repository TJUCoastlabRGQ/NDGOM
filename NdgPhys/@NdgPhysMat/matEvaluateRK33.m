function matEvaluateRK33( obj )

Nmesh = obj.Nmesh;
[rk0, rk1, rk2, rkt] = GetRKParamter();

time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys0 = cell( obj.Nmesh, 1 );
for n = 1:obj.Nmesh
    fphys0{n} = zeros( obj.meshUnion(n).cell.Np, obj.meshUnion(n).K, obj.Nvar );
end
fphys = obj.fphys;

visual = makeVisualizationFromNdgPhys( obj );
% init limiter and output file
hwait = waitbar(0,'Runing MatSolver....');
tic;
while( time < ftime )
    dt = 0.5 * obj.matUpdateTimeInterval( fphys );
%     display(dt);
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for n = 1:obj.Nmesh
        fphys0{n} = fphys{n};
    end
    
    for intRK = 1:3
        tloc = time + rkt(intRK) * dt;
        obj.matUpdateExternalField( tloc, fphys );
        obj.matEvaluateRHS( fphys );
        
        for n = 1:Nmesh
            fphys{n}(:,:, obj.varFieldIndex) ...
                = rk0(intRK)*fphys0{n}(:,:, obj.varFieldIndex)...
                + rk1(intRK)*fphys{n}(:,:, obj.varFieldIndex) ...
                + rk2(intRK)*dt*obj.frhs{n};
        end
        
%         data = [ 1, 1, 1 ];
%         fphys = obj.matEvaluateLimiter( fphys );
        fphys = obj.matEvaluatePostFunc( fphys );
        fphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata(obj, fphys, rk2(intRK)*dt);
    end
%     fprintf('processing %f...\n', time/ftime);
%     obj.draw( fphys );
    
    visual.drawResult( fphys{1}(:, :, 1) + fphys{1}(:, :, 4) ); 
    time = time + dt;
%     obj.matUpdateOutputResult( time, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ...
        ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
toc;
hwait.delete();
% obj.matUpdateFinalResult( time, fphys );
obj.fphys = fphys;
end

function [rk0, rk1, rk2, rkt] = GetRKParamter()
rk0 = [ 0.0, 3/4, 1/3 ];
rk1 = [ 1.0, 1/4, 2/3 ];
rk2 = [ 1.0, 1/4, 2/3 ];
rkt = [ 0.0, 1.0, 1.0];
end

