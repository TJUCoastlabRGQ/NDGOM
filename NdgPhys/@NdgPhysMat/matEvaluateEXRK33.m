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
%     dt = 0.2* obj.matUpdateTimeInterval( fphys );
    dt = 2.7;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    for n = 1:obj.Nmesh
        fphys0{n} = fphys{n};
    end
    
    tloc = time + ParameterC(1) * dt;
    obj.matUpdateExternalField( tloc, fphys );    
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
    fphys{1}(:,:,1) = fphys0{1}(:,:,1) + dt*EXb(1)*ExplicitRHS2d(:,:,1) +...
        dt*EXb(2)*ExplicitRHS2d(:,:,2) + dt*EXb(3)*ExplicitRHS2d(:,:,3) + ...
        dt*EXb(4)*ExplicitRHS2d(:,:,4);
    fphys{1}(:,:,2) = fphys0{1}(:,:,2) + dt*EXb(1)*ExplicitRHS2d(:,:,5) +...
        dt*EXb(2)*ExplicitRHS2d(:,:,6) + dt*EXb(3)*ExplicitRHS2d(:,:,7) + ...
        dt*EXb(4)*ExplicitRHS2d(:,:,8);
    fphys{1}(:,:,3) = fphys0{1}(:,:,3) + dt*EXb(1)*ExplicitRHS2d(:,:,9) +...
        dt*EXb(2)*ExplicitRHS2d(:,:,10) + dt*EXb(3)*ExplicitRHS2d(:,:,11) + ...
        dt*EXb(4)*ExplicitRHS2d(:,:,12);
    
    ExplicitRHS2d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K,4*obj.Nvar);
    
    visual.drawResult( fphys{1}(:, :, 1) ); 
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ...
        ['Runing MatSolver ', num2str( timeRatio ), '....']);    
end


end
function [Explicita, Explicitb, Parameterc] = GetRKParamter()
% data = roots([6 -18 9 -1]);
% GAMA = data(2);
% beta1 = -1.5*GAMA^2+4*GAMA-1/4;
% beta2 = 1.5*GAMA^2-5*GAMA+5/4;
% alpha1 = -0.35;
% alpha2 = (1/3-2*GAMA^2-2*beta2*alpha1*GAMA)/(GAMA*(1-GAMA));
GAMA = 0.435866521508459;
beta1 = 1.208496649176011;
beta2 = -0.644363170684470;
alpha1 = -0.35;
alpha2 = -0.989175724679849;
Explicita = [0 0 0 0;
    GAMA 0 0 0;
    (1+GAMA)/2-alpha1 alpha1 0 0;
    0 1-alpha2 alpha2 0];
Explicitb = [0 beta1 beta2 GAMA];
Parameterc = [0 GAMA (1+GAMA)/2 1];
end
