function matEvaluateEuler( obj )

Nmesh = obj.Nmesh;

time = obj.startTime;
% the external mode time
time2d = time;

ftime = obj.finalTime;

resQ2d = cell( obj.Nmesh, 1 );
resQ3d = cell( obj.Nmesh, 1 );
for n = 1:obj.Nmesh
    Np2 = obj.mesh2d(n).cell.Np;
    K2 = obj.mesh2d(n).K;
    resQ2d{n} = zeros( Np2, K2, obj.Nvar2d );
    
    Np3 = obj.mesh3d(n).cell.Np;
    K3 = obj.mesh3d(n).K;
    resQ3d{n} = zeros( Np3, K3, obj.Nvar );
end

fphys2d = obj.fphys2d;
fphys = obj.fphys;

visual = Visual2d( obj.mesh2d );

dt = obj.dt;
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    %     dt = obj.matUpdateTimeInterval( fphys2d );
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    while (time2d < time + dt)
        % time step for the two dimensional problem
        dt2d = 0.2 * obj.Solver2d.matUpdateTimeInterval( fphys2d );
        % check whether the two dimensional time step is larger than the three dimensional one or not
        if time2d + dt2d >= time + dt
            dt2d =  time + dt - time2d;
        end
        obj.Solver2d.matUpdateExternalField( time2d + dt2d, fphys2d );
        obj.Solver2d.matEvaluateRHS( fphys2d );
        obj.matEvaluate2dBoundaryStressRHS( obj.mesh3d, fphys, obj.Solver2d.frhs);
        for n = 1:Nmesh
            fphys2d{n}(:,:, obj.varFieldIndex2d) ...
                = fphys2d{n}(:,:, obj.varFieldIndex2d) + dt2d * obj.Solver2d.frhs{n};
        end
        
        %             fphys = obj.matEvaluatePostFunc( fphys );
        time2d = time2d + dt2d;
        visual.drawResult( fphys2d{1}(:, :, 1) );
    end
    %         tloc = time + rk4c( intRK ) * dt;
    %         obj.matUpdateExternalField( tloc, fphys2d, fphys3d );
    obj.matEvaluateRHS( fphys2d, fphys );
    
    for n = 1:Nmesh
        fphys{n}(:,:,1:2) = fphys{n}(:,:,1:2) + dt * obj.frhs{n};
        fphys{n} = obj.matEvaluateCorrectedMomentum( ...
            obj.mesh3d(n), fphys{n}, fphys2d{n} );
        fphys{n}(:,:,3) = obj.matEvaluateVerticalVelocity( ...
            obj.mesh3d(n), fphys2d{n}, fphys{n} );
    end
    time = time + dt;
    %     obj.matUpdateOutputResult( time, fphys2d, fphys );
    
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys2d = fphys2d;
obj.fphys = fphys;
obj.matUpdateFinalResult( time, fphys2d, fphys );
% obj.outputFile.closeOutputFile();
end