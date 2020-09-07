function matEvaluateSemiImplicit22(obj)
[rka, rkb, rkc, rkt] = GetRKParameter();
fphys2d = obj.fphys2d;
fphys = obj.fphys;
data = fphys{1}(:,:,1:2);
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
%The following three right hand side is not used actually
obj.ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,1);
obj.ExplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 2);
obj.ImplicitRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 2);
visual = Visual2d( obj.mesh2d );
hwait = waitbar(0,'Runing MatSolver....');
while( time < ftime )
    dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
    %       dt = 0.1;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,obj.varFieldIndex2d);
    Tempfphys = fphys{1}(:,:,obj.varFieldIndex);
    
    for intRK = 1:2
        tloc = time + rkt(intRK) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys2d, fphys );
        
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        fphys{1}(:,:,3) = obj.VertSolver.matCalculateVerticalVelocity( obj, fphys2d, fphys );
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
        obj.matCalculateExplicitRHSTerm( fphys2d, fphys, 1, 1);
        
        fphys2d{1}(:,:,1) = rka(intRK) * Tempfphys2d(:,:,1) + rkb(intRK) * fphys2d{1}(:,:,1) + rkc(intRK) * dt * obj.frhs2d{1}(:,:,1);
        fphys{1}(:,:,obj.varFieldIndex) = rka(intRK) * Tempfphys + rkb(intRK) * fphys{1}(:,:,obj.varFieldIndex) + rkc(intRK) * dt * obj.frhs{1}(:,:,1);
        
        %         [ fphys ] = obj.matImposeLimiter( fphys );
        %         [ fphys ] = obj.matFilterSolution( fphys );
        
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
        %         fphys2d = obj.matEvaluatePostFunc( fphys2d );
        
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    %dt here is problematic
    dt = 1;
%     [ fphys{1}(:,:,obj.varFieldIndex)] = ...
%         obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
%         10*ones(size(fphys2d{1}(:,:,1))), 10*ones(size(fphys{1}(:,:,4))), zeros(size(fphys{1}(:,:,obj.varFieldIndex))), dt );
    [ data ] = ...
        obj.VerticalEddyViscositySolver.matUpdateImplicitVerticalDiffusion( obj,...
        10*ones(size(fphys2d{1}(:,:,1))), 10*ones(size(fphys{1}(:,:,4))), data , dt );
    visual.drawResult( fphys2d{1}(:,:,1) );
    postplot(data(:,:,1));
    time = time + dt;
    obj.matUpdateOutputResult( time, fphys2d, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys2d = fphys2d;
obj.fphys = fphys;
obj.matUpdateFinalResult( time, fphys2d, fphys );
end

function [rka, rkb, rkc, rkt] = GetRKParameter
rka = [1 1/2];
rkb = [0 1/2];
rkc = [1 1/2];
rkt = [0 1];
end

function postplot(data)
figure;
index = [5,1];
Cor = [0,-0.1];
Np = 8;
Zcor = zeros(20,1);
Zcor(1:2)=Cor;
Dind = zeros(20,1);
Dind(1:2) = data(index);
for Layer = 2:10
    Zcor((Layer-1)*2+(1:2)) = Zcor((Layer-2)*2+(1:2)) - 0.1;
    Dind((Layer-1)*2+(1:2)) = data(index + (Layer-1)*Np);
end
plot(Dind, Zcor);
end