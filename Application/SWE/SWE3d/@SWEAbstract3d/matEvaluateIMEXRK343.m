function matEvaluateIMEXRK343( obj )
[EXa, IMa, EXb, IMb, ~] = GetRKParamter();
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');

fphys2d = obj.fphys2d;
fphys = obj.fphys;
%> allocate space for the rhs to be stored
ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,4);
ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 4*obj.Nvar);
ImplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
SystemRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nvar);
DiffusionCoefficient = fphys{1}(:,:,5)./fphys{1}(:,:,4).^2;
visual = Visual2d( obj.mesh2d );
dt = obj.dt;
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    %     dt = obj.matUpdateTimeInterval( fphys2d );
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,1);
    Tempfphys = fphys{1}(:,:,1:2);
    
    for intRK = 1:3
%         tloc = time + c( intRK ) * dt;
        %>Actually, boundary condition need to be imposed here
        %         obj.matUpdateExternalField( tloc, fphys2d, fphys );
        %This part need to consider the impact of the fext3d, as this is needed when impose the three-dimensional boundary
        [ExplicitRHS2d(:,:,intRK), ExplicitRHS3d(:,:,intRK), ExplicitRHS3d(:,:,intRK+4)] = ...
            matCalculateExplicitRHSTerm(obj, fphys2d, fphys, obj.fext2d);
        
        SystemRHS(:,:,1) = Tempfphys(:,:,1) + dt * EXa(intRK+1,1)*ExplicitRHS3d(:,:,1)+ dt * EXa(intRK+1,2)*ExplicitRHS3d(:,:,2)+...
            dt * EXa(intRK+1,3)*ExplicitRHS3d(:,:,3) + dt * IMa(intRK,1)*ImplicitRHS3d(:,:,1) + dt * IMa(intRK,2)*ImplicitRHS3d(:,:,2)+...
            dt * IMa(intRK,3)*ImplicitRHS3d(:,:,3);
        
        SystemRHS(:,:,2) = Tempfphys(:,:,2) + dt * EXa(intRK+1,1)*ExplicitRHS3d(:,:,5)+ dt * EXa(intRK+1,2)*ExplicitRHS3d(:,:,6)+...
            dt * EXa(intRK+1,3)*ExplicitRHS3d(:,:,7) + dt * IMa(intRK,1)*ImplicitRHS3d(:,:,4) + dt * IMa(intRK,2)*ImplicitRHS3d(:,:,5)+...
            dt * IMa(intRK,3)*ImplicitRHS3d(:,:,6);
        %Information about the 2d mesh is contained in meshUnion
        [ImplicitRHS3d(:,:,intRK), ImplicitRHS3d(:,:,intRK+3), fphys{1}(:,:,1), fphys{1}(:,:,2)] = ...
            matUpdateImplicitVerticalDiffusion(SystemRHS, DiffusionCoefficient, obj.meshUnion(1), fphys{1},IMa(intRK,intRK), dt, obj.Cf{1},...
            obj.WindTaux{1}, obj.WindTauy{1});
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );   
        %> update the vertical velocity
        fphys{1}(:,:,3) = obj.matEvaluateVerticalVelocity( obj.meshUnion(1), fphys2d{1}, fphys{1} );        
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
        % fphys2d = obj.matEvaluatePostFunc( fphys2d );
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    %>Actually, boundary condition need to be imposed here
    %     obj.matUpdateExternalField( time + dt, fphys2d, fphys );
    [ExplicitRHS2d(:,:,4), ExplicitRHS3d(:,:,4), ExplicitRHS3d(:,:,4+4)] = ...
        matCalculateExplicitRHSTerm(obj, fphys2d, fphys, obj.fext2d);
    %>Update the velocity
    fphys{1}(:,:,1) = Tempfphys(:,:,1) + dt * EXb(1) * ExplicitRHS3d(:,:,1) + dt * EXb(2) * ExplicitRHS3d(:,:,2)+...
        dt * EXb(3) * ExplicitRHS3d(:,:,3) + dt * EXb(4) * ExplicitRHS3d(:,:,4) + dt * IMb(1) * ImplicitRHS3d(:,:,1)+...
        dt * IMb(2) * ImplicitRHS3d(:,:,2) + dt * IMb(3) * ImplicitRHS3d(:,:,3);
    
    fphys{1}(:,:,2) = Tempfphys(:,:,2) + dt * EXb(1) * ExplicitRHS3d(:,:,5) + dt * EXb(2) * ExplicitRHS3d(:,:,6)+...
        dt * EXb(3) * ExplicitRHS3d(:,:,7) + dt * EXb(4) * ExplicitRHS3d(:,:,8) + dt * IMb(1) * ImplicitRHS3d(:,:,4)+...
        dt * IMb(2) * ImplicitRHS3d(:,:,5) + dt * IMb(3) * ImplicitRHS3d(:,:,6);
    
    fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXb(1) * ExplicitRHS2d(:,:,1) + dt * EXb(2) * ExplicitRHS2d(:,:,2)+...
        dt * EXb(3) * ExplicitRHS2d(:,:,3) + dt * EXb(4) * ExplicitRHS2d(:,:,4);
    visual.drawResult( fphys2d{1}(:,:,1) );
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
    ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,4);
    ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 4*obj.Nvar);
    ImplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
    time = time + dt;
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);       
%     [fphys{1}(:,:,5), obj.Cf{1}] = obj.EddyViscositySolver.matUpdateEddyViscosity( obj, obj.mesh2d, ...
%         obj.meshUnion(1), fphys2d, fphys, dt , time, obj.WindTaux{1}, obj.WindTauy{1} );
    
    %> Update the diffusion coefficient
    DiffusionCoefficient = fphys{1}(:,:,5)./fphys{1}(:,:,4).^2;
    display(time);
    obj.matUpdateOutputResult( time, fphys2d, fphys );
    timeRatio = time / ftime;
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys2d = fphys2d;
obj.fphys = fphys;
obj.matUpdateFinalResult( time, fphys2d, fphys );
% obj.outputFile.closeOutputFile();
end

function [ExplicitRHS2d, ExplicitHuRHS3d, ExplicitHvRHS3d] = matCalculateExplicitRHSTerm(obj, fphys2d, fphys, fext2d)
obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d, fphys, fext2d);
ExplicitRHS2d = obj.frhs2d{1}(:,:,1);
obj.advectionSolver.evaluateAdvectionRHS( fphys );
obj.matEvaluateSourceTerm( fphys );
ExplicitHuRHS3d = obj.frhs{1}(:,:,1);
ExplicitHvRHS3d = obj.frhs{1}(:,:,2);
end

function [ImplicithuRHS3d, ImplicithvRHS3d, hu, hv] = matUpdateImplicitVerticalDiffusion(SystemRHS, DiffusionCoefficient,...
    meshUnion, fphys, ImplicitParameter, dt, Cf, WindTaux, WindTauy)
%>Allocate memory space first
Tau = 30;
Nz = meshUnion.Nz;
ImplicithuRHS3d = zeros(size(meshUnion.x));
ImplicithvRHS3d = zeros(size(meshUnion.x));
hu              = zeros(size(meshUnion.x));
hv              = zeros(size(meshUnion.x));
StiffMatrix     = zeros(meshUnion.cell.Np*Nz);
BottomEidM      = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end-1)~=0,end-1);
UpEidM          = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end)~=0,end);
Np = meshUnion.cell.Np;
for i =1:meshUnion.mesh2d(1).K
    %> At present, we assume the mesh is uniform in the vertical direction
    ElementalMassMatrix3d = diag(meshUnion.J(:,(i-1)*Nz+1))*meshUnion.cell.M;
    ElementalMassMatrix2d = diag(meshUnion.mesh2d(1).J(:,i))*meshUnion.mesh2d(1).cell.M;
    Dz3d = diag(meshUnion.rz(:,(i-1)*Nz+1))*meshUnion.cell.Dr + diag(meshUnion.sz(:,(i-1)*Nz+1))*meshUnion.cell.Ds+...
        diag(meshUnion.tz(:,(i-1)*Nz+1))*meshUnion.cell.Dt;
    %> first cell first, then the other cells left
    LocalRows    = (1:Np)';
    AdjacentRows = (Np+1:2*Np)';
    LocalColumns = 1:Np;
    LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+1))*Dz3d;
    AdjacentPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+2))*Dz3d;
    %> Volume Integral Part
    OP11 = -Dz3d' * ElementalMassMatrix3d * LocalPhysicalDiffMatrix;
    %> Local Bottom Integral part
    OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP11);
    %> Adjacent bottom integral part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, AdjacentPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP12);
    %> Boundary part, not considered here
    [ SystemRHS(:,(i-1)*meshUnion.Nz + 1,1), SystemRHS(:,(i-1)*meshUnion.Nz + 1,2) ] = ImposeSurfaceBoundaryCondition(UpEidM, WindTaux(:,i),...
        WindTauy(:,i), ElementalMassMatrix2d, ElementalMassMatrix3d, dt, ImplicitParameter, SystemRHS(:,(i-1)*meshUnion.Nz + 1,1), SystemRHS(:,(i-1)*meshUnion.Nz + 1,2));
    %> Assemble into the StiffMatrix, it's noted that, the diagonal part has been included by eye(Np)
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
    StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
    for j = 2:meshUnion.Nz-1
        UpAdjacentRows = ((j-2)*Np+1:(j-1)*Np)';
        LocalRows    = ((j-1)*Np+1:j*Np)';
        BottomAdjacentRows = (j*Np+1:(j+1)*Np)';
        LocalColumns = (j-1)*Np+1:j*Np;
        UpPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+j-1))*Dz3d;
        LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+j))*Dz3d;
        BottomPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+j+1))*Dz3d;
        %> Volume Integral Part
        OP11 = -Dz3d' * ElementalMassMatrix3d * LocalPhysicalDiffMatrix;
        %> Local Bottom Integral part
        OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP11);
        %> Local Up Integral part
        OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP11);
        %> Assemble the local integral part into the StiffMatrix
        StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
        %> The upper adjacent cell part
        OP12 = zeros(meshUnion.cell.Np);
        OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, UpPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP12);
        StiffMatrix(UpAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
        %> The lower adjacent cell part
        OP12 = zeros(meshUnion.cell.Np);
        OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, BottomPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP12);
        StiffMatrix(BottomAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
    end
    %> for the bottom most cell
    AdjacentRows = ((meshUnion.Nz-2)*Np+1:(meshUnion.Nz-1)*Np)';
    LocalRows    = ((meshUnion.Nz-1)*Np+1:meshUnion.Nz*Np)';    
    LocalColumns = (meshUnion.Nz-1)*Np+1:meshUnion.Nz*Np;
    LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+Nz))*Dz3d;
    AdjacentPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+Nz-1))*Dz3d;
    %> Volume Integral Part
    OP11 = -Dz3d' * ElementalMassMatrix3d * LocalPhysicalDiffMatrix;
    %> Local Up Integral part
    OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP11);
    %> Impose bottom boundary condition
    OP11 = ImposeBottomBoundaryCondition(BottomEidM, OP11, ElementalMassMatrix2d, ...
        fphys(:,i*meshUnion.Nz,1), fphys(:,i*meshUnion.Nz,2), fphys(:,i*meshUnion.Nz,4), Cf);
    %> Assemble the local integral part into the StiffMatrix
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
    %> The upper adjacent cell part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, AdjacentPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP12);
    StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
    %This part is problematic, we need to consider the date structure
    temphuRHS = SystemRHS(:,(i-1)*meshUnion.Nz + 1:i*meshUnion.Nz,1);
    temphvRHS = SystemRHS(:,(i-1)*meshUnion.Nz + 1:i*meshUnion.Nz,2);
    hu((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np) = (eye(meshUnion.cell.Np*Nz) - ImplicitParameter * dt * StiffMatrix)\temphuRHS(:);
    hv((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np) = (eye(meshUnion.cell.Np*Nz) - ImplicitParameter * dt * StiffMatrix)\temphvRHS(:);
    %> This part is problematic, as the diagonal part has been included, this
    %> means that, the StiffMatrix is not the stiffmatrix assembled according to
    %> the weak form of the laplacian operator.
    ImplicithuRHS3d((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np) = StiffMatrix * hu((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np)';
    ImplicithvRHS3d((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np) = StiffMatrix * hv((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np)';
end
end

function OP11 = LocalUpBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
OP11(:, eidM)   = OP11(:, eidM)   + 0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d;
OP11(eidM, :)   = OP11(eidM, :)   + 0.5*massMatrix2d*Dz(eidM,:);
OP11(eidM,eidM) = OP11(eidM,eidM) - Tau*massMatrix2d;
end

function OP11 = LocalDownBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
OP11(:, eidM)   = OP11(:, eidM)   - 0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d;
OP11(eidM, :)   = OP11(eidM, :)   - 0.5*massMatrix2d*Dz(eidM,:);
OP11(eidM,eidM) = OP11(eidM,eidM) - Tau*massMatrix2d;
end

function OP12 = AdjacentDownBoundaryIntegral(eidM, eidP, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
%> Here, Down or up is relative to local cell
OP12(:,eidM)    = OP12(:,eidM) - 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;
OP12(eidP,:)    = OP12(eidP,:) + 0.5 * massMatrix2d * Dz(eidM,:);
OP12(eidP,eidM) = OP12(eidP,eidM) + Tau * massMatrix2d;
end

function OP12 = AdjacentUpBoundaryIntegral(eidM, eidP, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
OP12(:,eidM)    = OP12(:,eidM) + 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;
OP12(eidP,:)    = OP12(eidP,:) - 0.5 * massMatrix2d * Dz(eidM,:);
OP12(eidP,eidM) = OP12(eidP,eidM) + Tau * massMatrix2d;
end

function OP11 = ImposeBottomBoundaryCondition(eidM, OP11, massMatrix2d, hu, hv, h, Cf)
%> This part is positive, as this is teated implicitly
OP11(eidM, eidM) = OP11(eidM, eidM) + massMatrix2d * diag(Cf./h(eidM)./h(eidM).*sqrt(hu(eidM).^2 + hv(eidM).^2));
end

function [huRHS, hvRHS] = ImposeSurfaceBoundaryCondition(eidM, WindTaux, WindTauy, massMatrix2d, massMatrix3d, dt, ImplicitParameter, huRHS, hvRHS)
%> This part is negative, as this is teated explicitly
temphuRHS = zeros(size(WindTaux)); temphvRHS = zeros(size(WindTauy));
temphuRHS(eidM) = massMatrix2d*(dt*ImplicitParameter*WindTaux);
temphvRHS(eidM) = massMatrix2d*(dt*ImplicitParameter*WindTauy);
huRHS = huRHS - massMatrix3d\temphuRHS;
hvRHS = hvRHS - massMatrix3d\temphvRHS;
end

function [Explicita, Implicita, Explicitb, Implicitb, Parameterc] = GetRKParamter()
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
Implicita = [GAMA 0 0;
    (1-GAMA)/2 GAMA 0;
    beta1 beta2 GAMA];
Explicitb = [0 beta1 beta2 GAMA];
Implicitb = [beta1 beta2 GAMA];
end