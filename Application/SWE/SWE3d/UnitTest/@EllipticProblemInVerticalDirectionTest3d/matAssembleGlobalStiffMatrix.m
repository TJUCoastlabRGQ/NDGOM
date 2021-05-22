function [ StiffMatrix, LBRHS, RBRHS ] = matAssembleGlobalStiffMatrix(obj)
meshUnion = obj.meshUnion;
mesh2d = obj.meshUnion.mesh2d;
DiffusionCoefficient = ones(size(obj.meshUnion(1).x));
BottomEidM      = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end-1)~=0,end-1);
UpEidM          = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end)~=0,end);
%> Calculation of the penalty parameter
Tau = zeros(meshUnion.Nz+1, mesh2d.K);
Tau = matCalculatePenaltyParameter(meshUnion, mesh2d, DiffusionCoefficient, BottomEidM, UpEidM, Tau);
Nz = meshUnion.Nz;
Np = meshUnion.cell.Np;
StiffMatrix     = zeros(Np*Nz);
    %> At present, we assume the mesh is uniform in the vertical direction
ElementalMassMatrix3d = diag(meshUnion.J(:,1))*meshUnion.cell.M;
ElementalMassMatrix2d = diag(meshUnion.mesh2d(1).J(:,1))*meshUnion.mesh2d(1).cell.M;  
Dz3d = diag(meshUnion.rz(:,1))*meshUnion.cell.Dr + diag(meshUnion.sz(:,1))*meshUnion.cell.Ds+...
        diag(meshUnion.tz(:,1))*meshUnion.cell.Dt;
LocalRows    = (1:Np)';
AdjacentRows = (Np+1:2*Np)';
LocalColumns = 1:Np;
LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,1))*Dz3d;
AdjacentPhysicalDiffMatrix = diag(DiffusionCoefficient(:,2))*Dz3d;
OP11 = -Dz3d' * ElementalMassMatrix3d * LocalPhysicalDiffMatrix;
i = 1;
OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(2,i), OP11);
OP12 = zeros(meshUnion.cell.Np);
OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(2,i), OP12);
[ OP11, LBRHS ] = ImposeSurfaceDirichletBoundaryCondition(obj, UpEidM, LocalPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, ElementalMassMatrix3d, Tau(1,i), OP11);
StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
for j = 2:Nz-1
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
   [ RBRHS ] = ImposeBottomNewmannBoundaryCondition(obj, BottomEidM, ElementalMassMatrix2d, ElementalMassMatrix3d );
    %> Assemble the local integral part into the StiffMatrix
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
    %> The upper adjacent cell part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz3d, ElementalMassMatrix2d, Tau(Nz,i), OP12);
    StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
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
n0 = mesh2d.cell.Nface + 2;
%> here Nz stands for ratio between area of surface and volume of the studied cell
Nz = mesh3d.Nz;
for i = 1:mesh2d.K
    %> The surface most face for each column
    Tau(1,i) = (P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(UpEidM, (i-1)*Nz+1));
    for j = 2:Nz
        Tau(j,i) = (P+1)*(P+3)/3*n0/2*Nz*max(max(DiffusionCoefficient(BotEidM, (i-1)*Nz+j-1)),...
            max(DiffusionCoefficient(UpEidM, (i-1)*Nz+j)));
    end
    %> The bottom most face for each column
    Tau(Nz+1,i) = (P+1)*(P+3)/3*n0/2*Nz*max(DiffusionCoefficient(BotEidM, (i-1)*Nz+Nz));
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


function [ OP11, LBRHS ] = ImposeSurfaceDirichletBoundaryCondition(obj, UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, massMatrix3d, Tau, OP11)
% OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
epsilon = -1;
OP11 = DirichletTopBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
SurfTaux = obj.DirichExact(:,1);
temphuRHS = zeros(obj.meshUnion.cell.Np, 1);
%> move these terms about Dirichlet boundary to the left 
temphuRHS(UpEidM) = 1*Tau*massMatrix2d * SurfTaux;
temphuRHS = temphuRHS + epsilon * LocalPhysicalDiffMatrix(UpEidM,:)'*massMatrix2d*(-1*SurfTaux);
LBRHS = massMatrix3d\temphuRHS;
end

function OP11 = DirichletTopBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * 1 * physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   +  massMatrix2d*physicalDiffMatrix(eidM,:); %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d; %checked
end

function [ RBRHS ] = ImposeBottomNewmannBoundaryCondition(obj, eidM, massMatrix2d, massMatrix3d )
%> This part is negative, as this is teated explicitly
WindTaux = obj.NewmannExact;
temphuRHS = zeros(obj.meshUnion.cell.Np, 1);
temphuRHS(eidM) = massMatrix2d * WindTaux;
RBRHS = massMatrix3d\temphuRHS;
RBRHS = -1 * RBRHS;
end

