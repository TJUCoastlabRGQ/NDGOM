function matEvaluateEXRK33( obj )
[EXa, EXb, ParameterC] = GetRKParamter();

time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys0 = cell( obj.Nmesh, 1 );
for n = 1:obj.Nmesh
    fphys0{n} = zeros( obj.meshUnion(n).cell.Np, obj.meshUnion(n).K, obj.Nvar );
end
fphys = obj.fphys;

ExplicitRHS2d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K,4*obj.Nvar);

visual = makeVisualizationFromNdgPhys( obj );
% init limiter and output file
hwait = waitbar(0,'Runing MatSolver....');

while( time < ftime )
    dt = 0.2* obj.matUpdateTimeInterval( fphys );
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for n = 1:obj.Nmesh
        fphys0{n} = fphys{n};
    end
    
    obj.matEvaluateRHS( fphys );
    ExplicitRHS2d(:,:,1) = obj.frhs{1}(:,:,1);
    ExplicitRHS2d(:,:,1+4) = obj.frhs{1}(:,:,2);
    ExplicitRHS2d(:,:,1+4+4) = obj.frhs{1}(:,:,3);
    for intRK = 1:3
        tloc = time + ParameterC(intRK+1) * dt;
        obj.matUpdateExternalField( tloc, fphys );
        fphys{1}(:,:,1) = fphys0{1}(:,:,1) + dt * EXa(intRK+1,1) * ExplicitRHS2d(:,:,1) + ...
            dt * EXa(intRK+1,2) * ExplicitRHS2d(:,:,2)+ dt * EXa(intRK+1,3) * ExplicitRHS2d(:,:,3) + ...
            dt * EXa(intRK+1,4) * ExplicitRHS2d(:,:,4);
        fphys{1}(:,:,2) = fphys0{1}(:,:,2) + dt * EXa(intRK+1,1) * ExplicitRHS2d(:,:,5) + ...
            dt * EXa(intRK+1,2) * ExplicitRHS2d(:,:,6)+ dt * EXa(intRK+1,3) * ExplicitRHS2d(:,:,7) + ...
            dt * EXa(intRK+1,4) * ExplicitRHS2d(:,:,8);
        fphys{1}(:,:,3) = fphys0{1}(:,:,3) + dt * EXa(intRK+1,1) * ExplicitRHS2d(:,:,9) + ...
            dt * EXa(intRK+1,2) * ExplicitRHS2d(:,:,10)+ dt * EXa(intRK+1,3) * ExplicitRHS2d(:,:,11) + ...
            dt * EXa(intRK+1,4) * ExplicitRHS2d(:,:,12);
        obj.matEvaluateRHS( fphys );
        ExplicitRHS2d(:,:,intRK + 1) = obj.frhs{1}(:,:,1);
        ExplicitRHS2d(:,:,intRK + 1+4) = obj.frhs{1}(:,:,2);
        ExplicitRHS2d(:,:,intRK + 1+4+4) = obj.frhs{1}(:,:,3);
    end
    fphys{1}(:,:,1) = fphys0{1}(:,:,1) + EXb(1)*ExplicitRHS2d(:,:,1) +...
        EXb(2)*ExplicitRHS2d(:,:,2) + EXb(3)*ExplicitRHS2d(:,:,3) + ...
        EXb(4)*ExplicitRHS2d(:,:,4);
    fphys{1}(:,:,2) = fphys0{1}(:,:,2) + EXb(1)*ExplicitRHS2d(:,:,5) +...
        EXb(2)*ExplicitRHS2d(:,:,6) + EXb(3)*ExplicitRHS2d(:,:,7) + ...
        EXb(4)*ExplicitRHS2d(:,:,8);
    fphys{1}(:,:,3) = fphys0{1}(:,:,3) + EXb(1)*ExplicitRHS2d(:,:,9) +...
        EXb(2)*ExplicitRHS2d(:,:,10) + EXb(3)*ExplicitRHS2d(:,:,11) + ...
        EXb(4)*ExplicitRHS2d(:,:,12);
    
    visual.drawResult( fphys{1}(:, :, 1) ); 
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ...
        ['Runing MatSolver ', num2str( timeRatio ), '....']);    
end


end
function [Explicita, Explicitb, Parameterc] = GetRKParamter()
GAMA = 0.435866521508460;
beta1 = 1.208496649176012;
beta2 = -0.644363170684471;
alpha1 = -0.35;
alpha2 = -0.989175724679855;
Parameterc = [0 GAMA (1+GAMA)/2 1];
Explicita = [0 0 0 0;
    GAMA 0 0 0;
    (1+GAMA)/2-alpha1 alpha1 0 0;
    0 1-alpha2 alpha2 0];
Explicitb = [0 beta1 beta2 GAMA];
end
