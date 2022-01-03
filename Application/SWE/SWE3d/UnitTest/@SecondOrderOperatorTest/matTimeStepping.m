function matTimeStepping(obj)
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys2d = cell(1);
fphys2d{1} = zeros(size(obj.mesh2d(1).x));
fphys = obj.fphys;
%> allocate space for the rhs to be stored
ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 4*obj.Nvar);
ImplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
SystemRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K);
DiffusionCoefficient = obj.miu * ones(size(obj.meshUnion(1).x));
dt = 0.5;
% visual = Visual2d( obj.mesh2d );
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys = fphys{1}(:,:,1);
    
    for intRK = 1:3
        tloc = time + c( intRK + 1 ) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc );
        %> Calculate the right hand side for the global system about the three-dimensional horizontal momentum
        SystemRHS(:,:) = Tempfphys(:,:,1) + dt * EXa(intRK+1,1)*ExplicitRHS3d(:,:,1)+ dt * EXa(intRK+1,2)*ExplicitRHS3d(:,:,2)+...
            dt * EXa(intRK+1,3)*ExplicitRHS3d(:,:,3) + dt * IMa(intRK,1)*ImplicitRHS3d(:,:,1) + dt * IMa(intRK,2)*ImplicitRHS3d(:,:,2)+...
            dt * IMa(intRK,3)*ImplicitRHS3d(:,:,3);
        
        %> Calculating the right hand side corresponds to the discretization of the stiff diffusion term and return the intermediate horizontal momentum term
        [ImplicitRHS3d(:,:,intRK), fphys{1}(:,:,1)] = ...
            matUpdateImplicitVerticalDiffusion(obj, SystemRHS, DiffusionCoefficient, obj.meshUnion(1), obj.mesh2d(1), fphys{1},IMa(intRK,intRK), dt, tloc);
        
    end
    fphys{1}(:,:,1) = Tempfphys(:,:,1) + dt * EXb(1) * ExplicitRHS3d(:,:,1) + dt * EXb(2) * ExplicitRHS3d(:,:,2) + ...
        dt * EXb(3) * ExplicitRHS3d(:,:,3) + dt * EXb(4) * ExplicitRHS3d(:,:,4) + dt * IMb(1) * ImplicitRHS3d(:,:,1) + ...
        dt * IMb(2) * ImplicitRHS3d(:,:,2) + dt * IMb(3) * ImplicitRHS3d(:,:,3);
    
    
    ExplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 4*obj.Nvar);
    ImplicitRHS3d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3*obj.Nvar);
    time = time + dt;
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

function [ImplicithuRHS3d, hu] = matUpdateImplicitVerticalDiffusion(obj, SystemRHS, DiffusionCoefficient,...
    meshUnion, mesh2d, fphys, ImplicitParameter, dt, time)
BottomEidM      = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end-1)~=0,end-1);
UpEidM          = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end)~=0,end);
%> Calculation of the penalty parameter
Tau = zeros(meshUnion.Nz+1, mesh2d.K);
Tau = matCalculatePenaltyParameter(meshUnion, mesh2d, DiffusionCoefficient, BottomEidM, UpEidM, Tau);
Nz = meshUnion.Nz;
ImplicithuRHS3d = zeros(size(meshUnion.x));
% ImplicithvRHS3d = zeros(size(meshUnion.x));
hu              = zeros(size(meshUnion.x));
% hv              = zeros(size(meshUnion.x));
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
    [ SystemRHS(:,(i-1)*meshUnion.Nz + 1), SurfhuStiffMatrix ] = ImposeSurfaceNewmannBoundaryCondition(obj, UpEidM, time,...
        ElementalMassMatrix2d, ElementalMassMatrix3d, dt, ImplicitParameter, SystemRHS(:,(i-1)*meshUnion.Nz + 1));
    
    %> Impose surface Dirichlet boundary condition
    %     [ SystemRHS(:,(i-1)*meshUnion.Nz + 1), SurfhuStiffMatrix, OP11 ] = ImposeSurfaceDirichletBoundaryCondition(obj, UpEidM,...
    %         LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, ElementalMassMatrix3d, dt, ImplicitParameter, Tau(1,i), OP11, SystemRHS(:,(i-1)*meshUnion.Nz + 1), time);
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
    %> Impose bottom Dirichlet boundary condition
%     [ SystemRHS(:,i*meshUnion.Nz), BothuStiffMatirx, OP11 ] = ImposeBottomDirichletBoundaryCondition(obj, BottomEidM,...
%         LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, ElementalMassMatrix3d, dt, ImplicitParameter, Tau(Nz+1,i), OP11, SystemRHS(:,i*meshUnion.Nz), time);
    
    [ SystemRHS(:,i*meshUnion.Nz), BothuStiffMatirx ] = ImposeBottomNewmannBoundaryCondition(obj, BottomEidM, time,...
        ElementalMassMatrix2d, ElementalMassMatrix3d, dt, ImplicitParameter, SystemRHS(:,i*meshUnion.Nz));
    
    %> Assemble the local integral part into the StiffMatrix
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
    %> The upper adjacent cell part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(Nz,i), OP12);
    StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
    %This part is problematic, we need to consider the date structure
    temphuRHS = SystemRHS(:,(i-1)*meshUnion.Nz + 1:i*meshUnion.Nz,1);
    hu((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np) = (eye(meshUnion.cell.Np*Nz) - ImplicitParameter * dt * StiffMatrix)\temphuRHS(:);
    %> This part is problematic, as the diagonal part has been included, this
    %> means that, the StiffMatrix is not the stiffmatrix assembled according to
    %> the weak form of the laplacian operator.
    ImplicithuRHS3d((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np) = StiffMatrix * hu((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np)';
    %     ImplicithvRHS3d((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np) = StiffMatrix * hv((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np)';
    %> surface contribution to the discretization of the stiff term
    ImplicithuRHS3d((i-1)*meshUnion.Nz*Np+1:(i-1)*meshUnion.Nz*Np+Np) = ImplicithuRHS3d((i-1)*meshUnion.Nz*Np + 1 : (i-1)*meshUnion.Nz*Np+Np) + SurfhuStiffMatrix';
    ImplicithuRHS3d( i*meshUnion.Nz*Np-Np+1 : i*meshUnion.Nz*Np ) = ImplicithuRHS3d( i*meshUnion.Nz*Np-Np+1 : i*meshUnion.Nz*Np ) + BothuStiffMatirx';
    %> bottom contribution to the discretization of the stiff term
    %     ImplicithuRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) = ImplicithuRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) + BothuStiffMatrix';
    %> surface contribution to the discretization of the stiff term
    %     ImplicithvRHS3d((i-1)*meshUnion.Nz*Np+1:(i-1)*meshUnion.Nz*Np+Np) = ImplicithvRHS3d((i-1)*meshUnion.Nz*Np + 1 : (i-1)*meshUnion.Nz*Np+Np) + SurfhvStiffMatrix';
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
    Tau(1,i) = 10*(P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(UpEidM, (i-1)*Nz+1));
    for j = 2:Nz
        Tau(j,i) = 10*(P+1)*(P+3)/3*n0/2*Nz*max(max(DiffusionCoefficient(BotEidM, (i-1)*Nz+j-1)),...
            max(DiffusionCoefficient(UpEidM, (i-1)*Nz+j)));
    end
    %> The bottom most face for each column
    Tau(Nz+1,i) = 10*(P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(BotEidM, (i-1)*Nz+Nz));
end
end

function OP11 = LocalUpBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * 1 * 0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   +  0.5*massMatrix2d*physicalDiffMatrix(eidM,:); %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d; %checked
end

function OP11 = LocalDownBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * (-1) *  0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   -  0.5*massMatrix2d*physicalDiffMatrix(eidM,:);  %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d;   %checked
end

function OP12 = AdjacentDownBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
%> Here, Down or up is relative to local cell
epsilon = -1;
OP12(:,eidM)    = OP12(:,eidM) - epsilon * (-1) * 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;
OP12(eidP,:)    = OP12(eidP,:) +  0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);  %checked
OP12(eidP,eidM) = OP12(eidP,eidM) +  Tau * massMatrix2d;    %checked
end

function OP12 = AdjacentUpBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
epsilon = -1;
OP12(:,eidM)    = OP12(:,eidM) -  epsilon * (1) * 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;   %checked
OP12(eidP,:)    = OP12(eidP,:) -  0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);    %checked
OP12(eidP,eidM) = OP12(eidP,eidM) +  Tau * massMatrix2d;          %checked
end

function [huRHS,BothuStiffMatrix ,OP11 ] = ImposeBottomDirichletBoundaryCondition(obj, BottomEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, massMatrix3d, dt, ImplicitParameter, Tau, OP11, huRHS, time)
epsilon = -1;
% OP11 = LocalLeftBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
OP11 = DirichletBottomBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
BottomTaux = obj.DirichExact(:,2);
temphuRHS = zeros(size(huRHS));
temphuRHS(BottomEidM) = Tau*massMatrix2d * BottomTaux;
temphuRHS = temphuRHS + epsilon * LocalPhysicalDiffMatrix(BottomEidM,:)'*massMatrix2d*(-1*BottomTaux);
BothuStiffMatrix = massMatrix3d\temphuRHS;
huRHS = huRHS + dt * ImplicitParameter * BothuStiffMatrix;
end

function OP11 = DirichletBottomBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * (-1) * physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   -  massMatrix2d*physicalDiffMatrix(eidM,:);  %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d;   %checked
end

function [huRHS,SurfhuStiffMatrix ,OP11 ] = ImposeSurfaceDirichletBoundaryCondition(obj, UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, massMatrix3d, dt, ImplicitParameter, Tau, OP11, huRHS, time)
% OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
epsilon = -1;
OP11 = DirichletTopBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
SurfTaux = obj.DirichExact(:,1);
% miu = 0.01;
% SurfTaux = 1/sqrt(4*time+1)*exp(-(0+0.5).^2/miu/(4*time+1)) * ones(size(massMatrix2d,1),1);
temphuRHS = zeros(size(huRHS));
temphuRHS(UpEidM) = Tau*massMatrix2d * SurfTaux;
temphuRHS = temphuRHS + epsilon * LocalPhysicalDiffMatrix(UpEidM,:)'*massMatrix2d*(1*SurfTaux);
SurfhuStiffMatrix = massMatrix3d\temphuRHS;
huRHS = huRHS + dt * ImplicitParameter * SurfhuStiffMatrix;
end

function OP11 = DirichletTopBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * 1 * physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   +  massMatrix2d*physicalDiffMatrix(eidM,:); %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d; %checked
end

function [huRHS, huStiffMatrix] = ImposeSurfaceNewmannBoundaryCondition(obj, eidM, time, massMatrix2d, massMatrix3d, dt, ImplicitParameter, huRHS)

WindTaux = obj.NewmannExact(:,1);
temphuRHS = zeros(size(huRHS));
temphuRHS(eidM) = massMatrix2d * WindTaux*(1);
huStiffMatrix = massMatrix3d\temphuRHS;
% huRHS = huRHS - huStiffMatrix;
% hvRHS = hvRHS - hvStiffMatrix;
huRHS = huRHS + dt*ImplicitParameter*huStiffMatrix;
end

function [huRHS, huStiffMatrix] = ImposeBottomNewmannBoundaryCondition(obj, eidM, time, massMatrix2d, massMatrix3d, dt, ImplicitParameter, huRHS)

WindTaux = obj.NewmannExact(:,2);
temphuRHS = zeros(size(huRHS));
temphuRHS(eidM) = massMatrix2d * WindTaux*(-1);
huStiffMatrix = massMatrix3d\temphuRHS;
% huRHS = huRHS - huStiffMatrix;
% hvRHS = hvRHS - hvStiffMatrix;
huRHS = huRHS + dt*ImplicitParameter*huStiffMatrix;
end

function [Explicita, Implicita, Explicitb, Implicitb, Parameterc] = GetRKParamter()
data = roots([6 -18 9 -1]);
GAMA = data(2);
beta1 = -1.5*GAMA^2+4*GAMA-1/4;
beta2 = 1.5*GAMA^2-5*GAMA+5/4;
alpha1 = -0.35;
alpha2 = (1/3-2*GAMA^2-2*beta2*alpha1*GAMA)/(GAMA*(1-GAMA));
% GAMA = 0.435866521508460;
% beta1 = 1.208496649176012;
% beta2 = -0.644363170684471;
% alpha1 = -0.35;
% alpha2 = -0.989175724679855;
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