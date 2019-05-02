function matEvaluateHeun( obj )

Nmesh = obj.Nmesh;
[~, rk4b, rk4c] = GetRKParamter();

time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
resQ = cell( obj.Nmesh, 1 );

for n = 1:obj.Nmesh
    resQ{n} = zeros( obj.meshUnion(n).cell.Np, obj.meshUnion(n).K, obj.Nvar );
end
tempfphys = obj.fphys;
fphys = obj.fphys;

visual = makeVisualizationFromNdgPhys( obj );
% init limiter and output file
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = obj.matUpdateTimeInterval( fphys ) * 0.2;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for intRK = 1:2
        tloc = time + rk4c(intRK) * dt;
        obj.matUpdateExternalField( tloc, tempfphys );
        obj.matEvaluateRHS( tempfphys );
        
        for n = 1:Nmesh
            resQ{n} = dt*obj.frhs{n};
            tempfphys{n}(:,:, obj.varFieldIndex) ...
                = tempfphys{n}(:,:, obj.varFieldIndex) + rk4b(intRK)*resQ{n};
        end
        
%         tempfphys = obj.matEvaluateLimiter( tempfphys );
        tempfphys = obj.matEvaluatePostFunc( tempfphys );
        tempfphys = obj.NonhydrostaticSolver.NdgConservativeNonhydrostaticUpdata(obj, tempfphys,dt);
%         obj.meshUnion(1).draw( fphys{1}(:,:,1) );
%         drawnow;
    end
    
    fphys{1}(:,:,obj.varFieldIndex) = 1/2 * fphys{1}(:,:,obj.varFieldIndex) + ...
        1/2 * tempfphys{1}(:,:,obj.varFieldIndex);
    fphys{1}(:,:,6) = tempfphys{1}(:,:,6);
    tempfphys{1} = fphys{1};
    
    visual.drawResult( fphys{1}(:, :, 1) );
%     obj.meshUnion(1).draw( fphys{1}(:,:,1) );
%     drawnow;
    
%     fprintf('processing %f ...\n', time/ftime);
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ...
        ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.matUpdateFinalResult( time, fphys );
obj.fphys = fphys;
end

function [rk4a, rk4b, rk4c] = GetRKParamter()
rk4a = [ 0.0, 0.0];
rk4b = [ 1, 1 ];
rk4c = [ 1, 1 ];
end

