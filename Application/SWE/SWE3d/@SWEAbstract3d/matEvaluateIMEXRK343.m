function matEvaluateIMEXRK343( obj )
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
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
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    dt = 0.4 * obj.matUpdateTimeInterval( fphys2d );
%     dt = 2.7;
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys2d = fphys2d{1}(:,:,1);
    Tempfphys = fphys{1}(:,:,1:2);
    
    tloc = time + c( 1 ) * dt;
    obj.matUpdateExternalField( tloc, fphys2d, fphys );
    [ExplicitRHS2d(:,:,1), ExplicitRHS3d(:,:,1), ExplicitRHS3d(:,:,1+4)] = ...
        matCalculateExplicitRHSTerm(obj, fphys2d, fphys, obj.fext2d);
    
    for intRK = 1:3
        tloc = time + c( intRK + 1 ) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc, fphys2d, fphys );
        %> Calculate the intermediate water depth
        fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXa(intRK+1,1) * ExplicitRHS2d(:,:,1) + dt * EXa(intRK+1,2) * ExplicitRHS2d(:,:,2)+...
            dt * EXa(intRK+1,3) * ExplicitRHS2d(:,:,3) + dt * EXa(intRK+1,4) * ExplicitRHS2d(:,:,4);
        %> Calculate the right hand side for the global system about the three-dimensional horizontal momentum
        SystemRHS(:,:,1) = Tempfphys(:,:,1) + dt * EXa(intRK+1,1)*ExplicitRHS3d(:,:,1)+ dt * EXa(intRK+1,2)*ExplicitRHS3d(:,:,2)+...
            dt * EXa(intRK+1,3)*ExplicitRHS3d(:,:,3) + dt * IMa(intRK,1)*ImplicitRHS3d(:,:,1) + dt * IMa(intRK,2)*ImplicitRHS3d(:,:,2)+...
            dt * IMa(intRK,3)*ImplicitRHS3d(:,:,3);
        SystemRHS(:,:,2) = Tempfphys(:,:,2) + dt * EXa(intRK+1,1)*ExplicitRHS3d(:,:,5)+ dt * EXa(intRK+1,2)*ExplicitRHS3d(:,:,6)+...
            dt * EXa(intRK+1,3)*ExplicitRHS3d(:,:,7) + dt * IMa(intRK,1)*ImplicitRHS3d(:,:,4) + dt * IMa(intRK,2)*ImplicitRHS3d(:,:,5)+...
            dt * IMa(intRK,3)*ImplicitRHS3d(:,:,6);
        
        %> Calculating the right hand side corresponds to the discretization of the stiff diffusion term and return the intermediate horizontal momentum term 
        [ImplicitRHS3d(:,:,intRK), ImplicitRHS3d(:,:,intRK+3), fphys{1}(:,:,1), fphys{1}(:,:,2)] = ...
            matUpdateImplicitVerticalDiffusion(SystemRHS, DiffusionCoefficient, obj.meshUnion(1), obj.mesh2d(1), fphys{1},IMa(intRK,intRK), dt, obj.Cf{1},...
            obj.WindTaux{1}, obj.WindTauy{1});
        
        fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
        fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
        %> update the vertical velocity
        fphys{1}(:,:,3) = obj.matEvaluateVerticalVelocity( obj.meshUnion(1), fphys2d{1}, fphys{1} );
        
        fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
        fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
        %> Calculation of the right hand side corresponds to the discretization of the non-stiff term at stage intRK+1 with the
        %> newly calculated intermediate water depth and three-dimensional horizontal momentum
        [ExplicitRHS2d(:,:,intRK+1), ExplicitRHS3d(:,:,intRK+1), ExplicitRHS3d(:,:,intRK+1+4)] = ...
            matCalculateExplicitRHSTerm(obj, fphys2d, fphys, obj.fext2d);
        
        % fphys2d = obj.matEvaluateLimiter( fphys2d );
        fphys2d = obj.matEvaluatePostFunc( fphys2d );
        % visual.drawResult( fphys2d{1}(:,:,1) );
        % figure; obj.mesh3d.drawHorizonSlice( fphys3d{1}(:, :, 1) )
    end
    %> Update the water depth and the three-dimensional horizontal momentum at the next step
    fphys2d{1}(:,:,1) = Tempfphys2d(:,:,1) + dt * EXb(1) * ExplicitRHS2d(:,:,1) + dt * EXb(2) * ExplicitRHS2d(:,:,2) + ...
        dt * EXb(3) * ExplicitRHS2d(:,:,3) + dt * EXb(4) * ExplicitRHS2d(:,:,4);
    fphys{1}(:,:,1) = Tempfphys(:,:,1) + dt * EXb(1) * ExplicitRHS3d(:,:,1) + dt * EXb(2) * ExplicitRHS3d(:,:,2) + ...
        dt * EXb(3) * ExplicitRHS3d(:,:,3) + dt * EXb(4) * ExplicitRHS3d(:,:,4) + dt * IMb(1) * ImplicitRHS3d(:,:,1) + ...
        dt * IMb(2) * ImplicitRHS3d(:,:,2) + dt * IMb(3) * ImplicitRHS3d(:,:,3);
    fphys{1}(:,:,2) = Tempfphys(:,:,2) + dt * EXb(1) * ExplicitRHS3d(:,:,5) + dt * EXb(2) * ExplicitRHS3d(:,:,6) + ...
        dt * EXb(3) * ExplicitRHS3d(:,:,7) + dt * EXb(4) * ExplicitRHS3d(:,:,8) + dt * IMb(1) * ImplicitRHS3d(:,:,4) + ...
        dt * IMb(2) * ImplicitRHS3d(:,:,5) + dt * IMb(3) * ImplicitRHS3d(:,:,6);
    
    fphys2d{1}(:, :, 2) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 1) );
    fphys2d{1}(:, :, 3) = obj.meshUnion(1).VerticalColumnIntegralField( fphys{1}(:, :, 2) );
    %> update the vertical velocity
    fphys{1}(:,:,3) = obj.matEvaluateVerticalVelocity( obj.meshUnion(1), fphys2d{1}, fphys{1} );
    
    fphys{1}(: , :, 4) = obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) );
    fphys{1}(: , :, 7) = fphys{1}(: , :, 4) + fphys{1}(: , :, 6);
    
    %>Update the velocity
    visual.drawResult( fphys2d{1}(:,:,1) );
    % obj.drawVerticalSlice( 20, 1, fphys3d{1}(:, :, 3) * 1e7 );
    %> reallocate the space for the rhs
    ExplicitRHS2d = zeros(obj.mesh2d(1).cell.Np, obj.mesh2d(1).K,4);
    ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 4*obj.Nvar);
    ImplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
    time = time + dt;
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
    meshUnion, mesh2d, fphys, ImplicitParameter, dt, Cf, WindTaux, WindTauy)
BottomEidM      = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end-1)~=0,end-1);
UpEidM          = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end)~=0,end);
%> Calculation of the penalty parameter
Tau = zeros(meshUnion.Nz+1, mesh2d.K);
Tau = matCalculatePenaltyParameter(meshUnion, mesh2d, DiffusionCoefficient, BottomEidM, UpEidM, Tau);
Nz = meshUnion.Nz;
ImplicithuRHS3d = zeros(size(meshUnion.x));
ImplicithvRHS3d = zeros(size(meshUnion.x));
hu              = zeros(size(meshUnion.x));
hv              = zeros(size(meshUnion.x));
StiffMatrix     = zeros(meshUnion.cell.Np*Nz);
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
    OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(2,i), OP11);
    %> Adjacent bottom integral part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(2,i), OP12);
    %> Boundary part, not considered here
    [ SystemRHS(:,(i-1)*meshUnion.Nz + 1,1), SystemRHS(:,(i-1)*meshUnion.Nz + 1,2), SurfhuStiffMatrix, SurfhvStiffMatrix ] = ImposeSurfaceBoundaryCondition(UpEidM, WindTaux(:,i),...
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
        OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(j+1,i), OP11);
        %> Local Up Integral part
        OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(j,i), OP11);
        %> Assemble the local integral part into the StiffMatrix
        StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
        %> The upper adjacent cell part
        OP12 = zeros(meshUnion.cell.Np);
        OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(j,i), OP12);
        StiffMatrix(UpAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
        %> The lower adjacent cell part
        OP12 = zeros(meshUnion.cell.Np);
        OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(j+1,i), OP12);
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
    OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(Nz,i), OP11);
    %> Impose bottom boundary condition
    [ OP11 ] = ImposeBottomBoundaryCondition(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(Nz+1,i), OP11);
    %> Assemble the local integral part into the StiffMatrix
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
    %> The upper adjacent cell part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(Nz,i), OP12);
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
    %> surface contribution to the discretization of the stiff term
    ImplicithuRHS3d((i-1)*meshUnion.Nz*Np+1:(i-1)*meshUnion.Nz*Np+Np) = ImplicithuRHS3d((i-1)*meshUnion.Nz*Np + 1 : (i-1)*meshUnion.Nz*Np+Np) + SurfhuStiffMatrix';
    %> bottom contribution to the discretization of the stiff term
%     ImplicithuRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) = ImplicithuRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) + BothuStiffMatrix';
    %> surface contribution to the discretization of the stiff term
    ImplicithvRHS3d((i-1)*meshUnion.Nz*Np+1:(i-1)*meshUnion.Nz*Np+Np) = ImplicithvRHS3d((i-1)*meshUnion.Nz*Np + 1 : (i-1)*meshUnion.Nz*Np+Np) + SurfhvStiffMatrix'; 
    %> bottom contribution to the discretization of the stiff term
%     ImplicithvRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) = ImplicithvRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) + BothvStiffMatrix';        
end
end

function Tau = matCalculatePenaltyParameter(mesh3d, mesh2d, DiffusionCoefficient, BotEidM, UpEidM, Tau)
%> @brief Evaluating the penalty parameter used to penalize the jump between adjacet cell used in IPDG for second order operator
%>@detail In this version, the Interior Penalty Discontinuous Galerkin(IPDG) method is used to treat
%> the second order diffusion operator. To do so, the penalty parameter is calculated according to
%> [1] Shahbazi K. An explicit expression for the penalty parameter of the interior penalty method[J].
%> Journal of Computational Physics, 2005, 205(2): 401-407.
%> [2] Pestiaux A. Parameterization of subgrid-scale processes in finite element sea ice-ocean models[D].
%> UCL-Universit¨¦ Catholique de Louvain, 2015. pg:28.
%> The formula is '$\tau=\frac{(D_p+1)(D_p+d)}{d}\frac{n_0}{2}\frac{A}{V}\miu$'
%> @param[in] mesh3d The three-dimensional mesh object
%> @param[in] mesh2d The two-dimensional mesh object
%> @param[in] Tau The pre-allocated penalty parameter
%> @param[in] DiffusionCoefficient The scalar diffusion parameter
%> @param[out] Tau The calculated penalty parameter with size ( Nz + 1 ) * K2d
P = mesh2d.cell.N;
%> for prisms, number of faces is 5
n0 = 5;
%> here Nz stands for ratio between area of surface and volume of the studied cell
Nz = mesh3d.Nz;
for i = 1:mesh2d.K
    %> The surface most face for each column
%     Tau(1,i) = (P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(UpEidM, (i-1)*Nz+1));
    for j = 2:Nz
        Tau(j,i) = (P+1)*(P+3)/3*n0/2*Nz*max(max(DiffusionCoefficient(BotEidM, (i-1)*Nz+j-1)),...
            max(DiffusionCoefficient(UpEidM, (i-1)*Nz+j)));
    end
    %> The bottom most face for each column
    Tau(Nz+1,i) = (P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(BotEidM, (i-1)*Nz+Nz));
end
end

function OP11 = LocalUpBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
OP11(:, eidM)   = OP11(:, eidM)   + 0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d;
OP11(eidM, :)   = OP11(eidM, :)   + 0.5*massMatrix2d*physicalDiffMatrix(eidM,:);
OP11(eidM,eidM) = OP11(eidM,eidM) - Tau*massMatrix2d;
end

function OP11 = LocalDownBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
OP11(:, eidM)   = OP11(:, eidM)   - 0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d;
OP11(eidM, :)   = OP11(eidM, :)   - 0.5*massMatrix2d*physicalDiffMatrix(eidM,:);
OP11(eidM,eidM) = OP11(eidM,eidM) - Tau*massMatrix2d;
end

function OP12 = AdjacentDownBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
%> Here, Down or up is relative to local cell
OP12(:,eidM)    = OP12(:,eidM) - 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;
OP12(eidP,:)    = OP12(eidP,:) + 0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);
OP12(eidP,eidM) = OP12(eidP,eidM) + Tau * massMatrix2d;
end

function OP12 = AdjacentUpBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
OP12(:,eidM)    = OP12(:,eidM) + 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;
OP12(eidP,:)    = OP12(eidP,:) - 0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);
OP12(eidP,eidM) = OP12(eidP,eidM) + Tau * massMatrix2d;
end

function [ OP11 ] = ImposeBottomBoundaryCondition(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP11)
OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau, OP11);
end

function [huRHS, hvRHS, huStiffMatrix, hvStiffMatrix] = ImposeSurfaceBoundaryCondition(eidM, WindTaux, WindTauy, massMatrix2d, massMatrix3d, dt, ImplicitParameter, huRHS, hvRHS)
%> This part is negative, as this is teated explicitly
temphuRHS = zeros(size(WindTaux)); temphvRHS = zeros(size(WindTauy));
temphuRHS(eidM) = massMatrix2d*(dt*ImplicitParameter*WindTaux);
temphvRHS(eidM) = massMatrix2d*(dt*ImplicitParameter*WindTauy);
huStiffMatrix = massMatrix3d\temphuRHS;
hvStiffMatrix = massMatrix3d\temphvRHS;
huRHS = huRHS - huStiffMatrix;
hvRHS = hvRHS - hvStiffMatrix;
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