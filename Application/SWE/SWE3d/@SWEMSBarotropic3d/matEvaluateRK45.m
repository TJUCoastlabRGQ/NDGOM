function matEvaluateRK45( obj )

Nmesh = obj.Nmesh;
[rk4a, rk4b, rk4c] = GetRKParamter();

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
        dt2d = obj.Solver2d.matUpdateTimeInterval( fphys2d );
        % check whether the two dimensional time step is larger than the three dimensional one or not
        if time2d + dt2d >= time + dt
            dt2d =  time + dt - time2d;
        end
        for intRK = 1:5
            tloc = time2d + rk4c( intRK ) * dt2d;
            obj.Solver2d.matUpdateExternalField( tloc, fphys2d );
            obj.Solver2d.matEvaluateRHS( fphys2d );
            obj.matEvaluate2dBoundaryStressRHS( obj.mesh3d, fphys, obj.Solver2d.frhs);
            for n = 1:Nmesh
                resQ2d{n} = rk4a( intRK ) * resQ2d{n} + dt2d * obj.Solver2d.frhs{n};
                fphys2d{n}(:,:, obj.varFieldIndex2d) ...
                    = fphys2d{n}(:,:, obj.varFieldIndex2d) + rk4b(intRK) * resQ2d{n};
            end
            
%             fphys = obj.matEvaluatePostFunc( fphys );
        end
        time2d = time2d + dt2d;
        visual.drawResult( fphys2d{1}(:, :, 1) );        
    end
    for intRK = 1:5
%         tloc = time + rk4c( intRK ) * dt;
        %         obj.matUpdateExternalField( tloc, fphys2d, fphys3d );
        obj.matEvaluateRHS( fphys2d, fphys );
        
        for n = 1:Nmesh
            resQ3d{n} = rk4a( intRK ) * resQ3d{n} + dt * obj.frhs{n};
            fphys{n}(:,:,1:2) = fphys{n}(:,:,1:2) + rk4b(intRK) * resQ3d{n};
            fphys{n} = obj.matEvaluateCorrectedMomentum( ...
                obj.mesh3d(n), fphys{n}, fphys2d{n} );
            fphys{n}(:,:,3) = obj.matEvaluateVerticalVelocity( ...
                obj.mesh3d(n), fphys2d{n}, fphys{n} );
        end
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

function [rk4a, rk4b, rk4c] = GetRKParamter()

rk4a = [            0.0 ...
    -567301805773.0/1357537059087.0 ...
    -2404267990393.0/2016746695238.0 ...
    -3550918686646.0/2091501179385.0  ...
    -1275806237668.0/842570457699.0];
rk4b = [ 1432997174477.0/9575080441755.0 ...
    5161836677717.0/13612068292357.0 ...
    1720146321549.0/2090206949498.0  ...
    3134564353537.0/4481467310338.0  ...
    2277821191437.0/14882151754819.0];
rk4c = [             0.0  ...
    1432997174477.0/9575080441755.0 ...
    2526269341429.0/6820363962896.0 ...
    2006345519317.0/3224310063776.0 ...
    2802321613138.0/2924317926251.0];

end