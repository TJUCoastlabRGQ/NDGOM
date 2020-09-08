function matTimeStepping222(obj)
[EXa, IMa, EXb, IMb, c] = GetRKParamter();
time = obj.getOption('startTime');
ftime = obj.getOption('finalTime');
fphys = obj.fphys;
%> allocate space for the rhs to be stored
ExplicitRHS1d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3);
ImplicitRHS1d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 2);
SystemRHS = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K);
DiffusionCoefficient = obj.miu * ones(size(obj.meshUnion(1).x));
dt = 0.005;
% visual = Visual2d( obj.mesh2d );
visual = makeVisualizationFromNdgPhys( obj );
hwait = waitbar(0,'Runing MatSolver....');
% try
while( time < ftime )
    if( time + dt > ftime )
        dt = ftime - time;
    end
    
    Tempfphys = fphys{1}(:,:,1);
    
    for intRK = 1:2
        tloc = time + c( intRK+1 ) * dt;
        %>Actually, boundary condition need to be imposed here
        obj.matUpdateExternalField( tloc );        
        %> Calculate the right hand side for the global system about the three-dimensional horizontal momentum
        SystemRHS(:,:) = Tempfphys + dt * EXa(intRK+1,1)*ExplicitRHS1d(:,:,1)+ dt * EXa(intRK+1,2)*ExplicitRHS1d(:,:,2)+...
            dt * EXa(intRK+1,3)*ExplicitRHS1d(:,:,3) + dt * IMa(intRK,1)*ImplicitRHS1d(:,:,1) + dt * IMa(intRK,2)*ImplicitRHS1d(:,:,2);
        
        %> Calculating the right hand side corresponds to the discretization of the stiff diffusion term and return the intermediate horizontal momentum term 
        [ImplicitRHS1d(:,:,intRK), fphys{1}(:,:,1)] = ...
            matUpdateImplicitVerticalDiffusion(obj, SystemRHS, DiffusionCoefficient, obj.meshUnion(1), fphys{1},IMa(intRK,intRK), dt, tloc);
        
    end
    fphys{1}(:,:,1) = Tempfphys(:,:,1) + dt * EXb(1) * ExplicitRHS1d(:,:,1) + dt * EXb(2) * ExplicitRHS1d(:,:,2) + ...
        dt * EXb(3) * ExplicitRHS1d(:,:,3) + dt * IMb(1) * ImplicitRHS1d(:,:,1) + ...
        dt * IMb(2) * ImplicitRHS1d(:,:,2);
    
    
    ExplicitRHS1d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 3);
    ImplicitRHS1d = zeros(obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, 2);
    time = time + dt;
    display(time);
%     obj.matUpdateOutputResult( time, fphys );
    timeRatio = time / ftime;
    visual.drawResult( fphys{1}(:, :, 1) );
    waitbar( timeRatio, hwait, ['Runing MatSolver ', num2str( timeRatio ), '....']);
end
hwait.delete();
obj.fphys = fphys;
% obj.matUpdateFinalResult( time, fphys );
% obj.outputFile.closeOutputFile();
end

function [ImplicithuRHS3d, hu] = matUpdateImplicitVerticalDiffusion(obj, SystemRHS, DiffusionCoefficient,...
    meshUnion, fphys, ImplicitParameter, dt, time)
LeftEidM      = meshUnion.cell.Fmask(1);
RightEidM          = meshUnion.cell.Fmask(2);
%> Calculation of the penalty parameter
Tau = zeros(1, meshUnion.K);
Tau = matCalculatePenaltyParameter(meshUnion, DiffusionCoefficient, LeftEidM, RightEidM, Tau);
Nz = meshUnion.K;
ImplicithuRHS3d = zeros(size(meshUnion.x));
% ImplicithvRHS3d = zeros(size(meshUnion.x));
hu              = zeros(size(meshUnion.x));
% hv              = zeros(size(meshUnion.x));
StiffMatrix     = zeros(meshUnion.cell.Np*Nz);
Np = meshUnion.cell.Np;
    %> At present, we assume the mesh is uniform in the vertical direction
    ElementalMassMatrix = diag(meshUnion.J(:,1))*meshUnion.cell.M;
    ElementalMassMatrix2d = 1;
    Dz1d = diag(meshUnion.rx(:,1))*meshUnion.cell.Dr;
    %> first cell first, then the other cells left
    LocalRows    = (1:Np)';
    AdjacentRows = (Np+1:2*Np)';
    LocalColumns = 1:Np;
    LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,1)) * Dz1d;
    AdjacentPhysicalDiffMatrix = diag(DiffusionCoefficient(:,2)) * Dz1d;
    %> Volume Integral Part
    OP11 = -Dz1d' * ElementalMassMatrix * LocalPhysicalDiffMatrix;
    %> Local Bottom Integral part
    OP11 = LocalRightBoundaryIntegral(RightEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(2), OP11);
    %> Adjacent bottom integral part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentRightBoundaryIntegral(RightEidM, LeftEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(2), OP12);
    %> Boundary part, not considered here
    [ SystemRHS(:,1), SurfhuStiffMatrix ] = ImposeSurfaceNewmannBoundaryCondition(obj, LeftEidM, time,...
        ElementalMassMatrix2d, ElementalMassMatrix, dt, ImplicitParameter, SystemRHS(:,1));
%     [ SystemRHS(:,1), SurfhuStiffMatrix, OP11 ] = ImposeSurfaceDirichletBoundaryCondition(obj, LeftEidM,...
%         LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, ElementalMassMatrix, dt, ImplicitParameter, Tau(1), OP11, SystemRHS(:,1), time);
    
    %> Assemble into the StiffMatrix, it's noted that, the diagonal part has been included by eye(Np)
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix\OP11;
    StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
    for j = 2:meshUnion.K-1
        UpAdjacentRows = ((j-2)*Np+1:(j-1)*Np)';
        LocalRows    = ((j-1)*Np+1:j*Np)';
        BottomAdjacentRows = (j*Np+1:(j+1)*Np)';
        LocalColumns = (j-1)*Np+1:j*Np;
        UpPhysicalDiffMatrix = diag(DiffusionCoefficient(:,j-1)) * Dz1d;
        LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,j)) * Dz1d;
        BottomPhysicalDiffMatrix = diag(DiffusionCoefficient(:,j+1)) * Dz1d;
        %> Volume Integral Part
        OP11 = -Dz1d' * ElementalMassMatrix * LocalPhysicalDiffMatrix;
        %> Local Bottom Integral part
        OP11 = LocalLeftBoundaryIntegral(LeftEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP11);
        %> Local Up Integral part
        OP11 = LocalRightBoundaryIntegral(RightEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP11);
        %> Assemble the local integral part into the StiffMatrix
        StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix\OP11;
        %> The upper adjacent cell part
        OP12 = zeros(meshUnion.cell.Np);
        OP12 = AdjacentRightBoundaryIntegral(RightEidM, LeftEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP12);
        StiffMatrix(BottomAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
        %> The lower adjacent cell part
        OP12 = zeros(meshUnion.cell.Np);
        OP12 = AdjacentLeftBoundaryIntegral(LeftEidM, RightEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP12);
        StiffMatrix(UpAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
    end
    %> for the bottom most cell
    AdjacentRows = ((Nz-2)*Np+1:(Nz-1)*Np)';
    LocalRows    = ((Nz-1)*Np+1:Nz*Np)';
    LocalColumns = (Nz-1)*Np+1:Nz*Np;
    LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,Nz)) * Dz1d;
    AdjacentPhysicalDiffMatrix = diag(DiffusionCoefficient(:,Nz-1)) * Dz1d;
    %> Volume Integral Part
    OP11 = -Dz1d' * ElementalMassMatrix * LocalPhysicalDiffMatrix;
    %> Local Up Integral part
    OP11 = LocalLeftBoundaryIntegral(LeftEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP11);
    %> Impose bottom Dirichlet boundary condition
    [ SystemRHS(:,Nz), BothuStiffMatirx, OP11 ] = ImposeBottomDirichletBoundaryCondition(obj, RightEidM,...
        LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, ElementalMassMatrix, dt, ImplicitParameter, Tau(1), OP11, SystemRHS(:,Nz), time);
%         [ SystemRHS(:,Nz), BothuStiffMatirx ] = ImposeBottomNewmannBoundaryCondition(obj, RightEidM, time,...
%         ElementalMassMatrix2d, ElementalMassMatrix, dt, ImplicitParameter, SystemRHS(:,Nz));
    %> Assemble the local integral part into the StiffMatrix
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix\OP11;
    %> The upper adjacent cell part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentLeftBoundaryIntegral(LeftEidM, RightEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP12);
    StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
    %This part is problematic, we need to consider the date structure
    temphuRHS = SystemRHS(:,:,1);
    hu(1:Nz*Np) = (eye(Np*Nz) - ImplicitParameter * dt * StiffMatrix)\temphuRHS(:);
    %> This part is problematic, as the diagonal part has been included, this
    %> means that, the StiffMatrix is not the stiffmatrix assembled according to
    %> the weak form of the laplacian operator.
    ImplicithuRHS3d(1:Nz*Np) = StiffMatrix * hu(1:Nz*Np)';
%     ImplicithvRHS3d((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np) = StiffMatrix * hv((i-1)*meshUnion.Nz*Np + 1 : i*meshUnion.Nz*Np)';
    %> surface contribution to the discretization of the stiff term
    ImplicithuRHS3d(1:Np) = ImplicithuRHS3d(1:Np) + SurfhuStiffMatrix';
    ImplicithuRHS3d( Nz*Np-Np+1:Nz*Np ) = ImplicithuRHS3d(Nz*Np-Np+1:Nz*Np ) + BothuStiffMatirx';
    %> bottom contribution to the discretization of the stiff term
%     ImplicithuRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) = ImplicithuRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) + BothuStiffMatrix';
    %> surface contribution to the discretization of the stiff term
%     ImplicithvRHS3d((i-1)*meshUnion.Nz*Np+1:(i-1)*meshUnion.Nz*Np+Np) = ImplicithvRHS3d((i-1)*meshUnion.Nz*Np + 1 : (i-1)*meshUnion.Nz*Np+Np) + SurfhvStiffMatrix'; 
    %> bottom contribution to the discretization of the stiff term
%     ImplicithvRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) = ImplicithvRHS3d((i-1)*meshUnion.Nz*Np+(meshUnion.Nz-1)*Np+1:i*meshUnion.Nz*Np) + BothvStiffMatrix';        
end

function Tau = matCalculatePenaltyParameter(mesh, DiffusionCoefficient, LeftEidM, RightEidM, Tau)
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
P = mesh.cell.N;
%> for prisms, number of faces is 5
n0 = 2;
%> here Nz stands for ratio between area of surface and volume of the studied cell
% Nz = mesh3d.Nz;
for i = 1:numel(Tau)
    Tau(i) = (P+1)*(P+1)/1*n0/2*1/mesh.LAV(1)*DiffusionCoefficient(1);
end
% for i = 1:mesh2d.K
%     %> The surface most face for each column
%     Tau(1,i) = 10*(P+1)*(P+1)/1*n0/2*Nz*max(DiffusionCoefficient(LeftEidM, (i-1)*Nz+1));
%     for j = 2:Nz
%         Tau(j,i) = 10*(P+1)*(P+3)/3*n0/2*Nz*max(max(DiffusionCoefficient(BotEidM, (i-1)*Nz+j-1)),...
%             max(DiffusionCoefficient(UpEidM, (i-1)*Nz+j)));
%     end
%     %> The bottom most face for each column
%     Tau(Nz+1,i) = 10*(P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(BotEidM, (i-1)*Nz+Nz));
% end
end

function OP11 = LocalRightBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * 1 * 0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   +  0.5*massMatrix2d*physicalDiffMatrix(eidM,:); %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d; %checked
end

function OP11 = LocalLeftBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * (-1) *  0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   -  0.5*massMatrix2d*physicalDiffMatrix(eidM,:);  %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d;   %checked
end

function OP12 = AdjacentLeftBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
%> Here, Down or up is relative to local cell
epsilon = -1;
OP12(:,eidM)    = OP12(:,eidM) - epsilon * (-1) * 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;   
OP12(eidP,:)    = OP12(eidP,:) +  0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);  %checked
OP12(eidP,eidM) = OP12(eidP,eidM) +  Tau * massMatrix2d;    %checked
end

function OP12 = AdjacentRightBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
epsilon = -1;
OP12(:,eidM)    = OP12(:,eidM) -  epsilon * (1) * 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;   %checked
OP12(eidP,:)    = OP12(eidP,:) -  0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);    %checked
OP12(eidP,eidM) = OP12(eidP,eidM) +  Tau * massMatrix2d;          %checked
end

function [huRHS,BothuStiffMatrix ,OP11 ] = ImposeBottomDirichletBoundaryCondition(obj, RightEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, massMatrix3d, dt, ImplicitParameter, Tau, OP11, huRHS, time)
epsilon = -1;
% OP11 = LocalRightBoundaryIntegral(RightEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
OP11 = DirichletRightBoundaryCondition(RightEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
BottomTaux = obj.DirichExact(2);
temphuRHS = zeros(size(huRHS));
temphuRHS(RightEidM) = Tau*massMatrix2d * BottomTaux;
temphuRHS = temphuRHS + epsilon * LocalPhysicalDiffMatrix(RightEidM,:)'*massMatrix2d*(1*BottomTaux);
BothuStiffMatrix = massMatrix3d\temphuRHS;
huRHS = huRHS + dt * ImplicitParameter * BothuStiffMatrix;
end

function OP11 = DirichletRightBoundaryCondition(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * 1 * physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   +  massMatrix2d*physicalDiffMatrix(eidM,:); %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d; %checked
end

function [huRHS,SurfhuStiffMatrix ,OP11 ] = ImposeSurfaceDirichletBoundaryCondition(obj, UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, massMatrix3d, dt, ImplicitParameter, Tau, OP11, huRHS, time)
epsilon = -1;
% OP11 = LocalLeftBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
OP11 = DirichletLeftBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);

SurfTaux = obj.DirichExact(1);
temphuRHS = zeros(size(huRHS));
temphuRHS(UpEidM) = Tau*massMatrix2d * SurfTaux;
temphuRHS = temphuRHS + epsilon * LocalPhysicalDiffMatrix(UpEidM,:)'*massMatrix2d*(-1*SurfTaux);
SurfhuStiffMatrix = massMatrix3d\temphuRHS;
huRHS = huRHS + dt * ImplicitParameter * SurfhuStiffMatrix;
end

function OP11 = DirichletLeftBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * (-1) * physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   -  massMatrix2d*physicalDiffMatrix(eidM,:);  %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d;   %checked
end

function [huRHS, huStiffMatrix] = ImposeBottomNewmannBoundaryCondition(obj, eidM, time, massMatrix2d, massMatrix, dt, ImplicitParameter, huRHS)
%> This part is negative, as this is teated explicitly
WindTaux = obj.NewmannExact(2);
temphuRHS = zeros(size(huRHS)); 
temphuRHS(eidM) = massMatrix2d * WindTaux*(1);
huStiffMatrix = massMatrix\temphuRHS;
% huRHS = huRHS - huStiffMatrix;
% hvRHS = hvRHS - hvStiffMatrix;
huRHS = huRHS + dt*ImplicitParameter*huStiffMatrix;
end

function [huRHS, huStiffMatrix] = ImposeSurfaceNewmannBoundaryCondition(obj, eidM, time, massMatrix2d, massMatrix, dt, ImplicitParameter, huRHS)
%> This part is negative, as this is teated explicitly
WindTaux = obj.NewmannExact(1);
temphuRHS = zeros(size(huRHS)); 
temphuRHS(eidM) = massMatrix2d * WindTaux*(-1);
huStiffMatrix = massMatrix\temphuRHS;
% huRHS = huRHS - huStiffMatrix;
% hvRHS = hvRHS - hvStiffMatrix;
huRHS = huRHS + dt*ImplicitParameter*huStiffMatrix;
end



function [Explicita, Implicita, Explicitb, Implicitb, Parameterc] = GetRKParamter()
GAMA = (2-sqrt(2))/2;
delta = 1-1/(2*GAMA);
Parameterc = [0 GAMA 1];
Explicita = [0 0 0;
    GAMA 0 0;
    delta 1-delta 0];
Implicita = [GAMA 0;
    (1-GAMA) GAMA];
Explicitb = [delta 1-delta 0];
Implicitb = [1-GAMA GAMA];
end