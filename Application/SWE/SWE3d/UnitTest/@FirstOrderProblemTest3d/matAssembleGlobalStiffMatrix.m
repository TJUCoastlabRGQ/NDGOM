function [ StiffMatrix, LBRHS ] = matAssembleGlobalStiffMatrix(obj)

% StiffMatrix = obj.NonhydrostaticSolver.PNPS;
%
% LBRHS = zeros(obj.mesh3d.cell.Np, 1);

meshUnion = obj.meshUnion;
BottomEidM      = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end-1)~=0,end-1);
UpEidM          = meshUnion.cell.Fmask(meshUnion.cell.Fmask(:,end)~=0,end);
%> Calculation of the penalty parameter
K = meshUnion.K;
Nz = meshUnion.Nz;
Np = meshUnion.cell.Np;
StiffMatrix     = zeros(K*Np);
ElementalMassMatrix3d = diag(meshUnion.J(:,1))*meshUnion.cell.M;
ElementalMassMatrix2d = diag(meshUnion.mesh2d(1).J(:,1))*meshUnion.mesh2d(1).cell.M;
Dz3d = diag(meshUnion.rz(:,1))*meshUnion.cell.Dr + diag(meshUnion.sz(:,1))*meshUnion.cell.Ds+...
    diag(meshUnion.tz(:,1))*meshUnion.cell.Dt;
for k = 1:meshUnion.mesh2d(1).K
    LocalRows    = ((k-1)*Np*Nz + 1:(k-1)*Np*Nz + Np)';
    AdjacentRows = ((k-1)*Np*Nz + Np+1:(k-1)*Np*Nz + 2*Np)';
    LocalColumns = (k-1)*Np*Nz + 1:(k-1)*Np*Nz + Np;
    
    OP11 = ElementalMassMatrix3d * Dz3d;
    
    OP11 = LocalDownBoundaryIntegral(BottomEidM, ElementalMassMatrix2d, OP11);
    
    OP11 = LocalDownLDGBoundaryIntegral(BottomEidM, ElementalMassMatrix2d, OP11);
    
    [ OP11, LBRHS ] = ImposeSurfaceDirichletBoundaryCondition(obj, UpEidM, ElementalMassMatrix3d, ElementalMassMatrix2d, OP11);
    
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, ElementalMassMatrix2d, OP12);
    
    OP12 = AdjacentDownLDGBoundaryIntegral(BottomEidM, UpEidM, ElementalMassMatrix2d, OP12);
    
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
    StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
    for j = 2:Nz-1
        UpAdjacentRows = ((k-1)*Np*Nz + (j-2)*Np+1:(k-1)*Np*Nz + (j-1)*Np)';
        LocalRows    = ((k-1)*Np*Nz + (j-1)*Np+1:(k-1)*Np*Nz + j*Np)';
        BottomAdjacentRows = ((k-1)*Np*Nz + j*Np+1:(k-1)*Np*Nz + (j+1)*Np)';
        LocalColumns = (k-1)*Np*Nz + (j-1)*Np+1:(k-1)*Np*Nz + j*Np;
        %> Volume Integral Part
        OP11 = ElementalMassMatrix3d * Dz3d;
        %> Local Bottom Integral part
        OP11 = LocalDownBoundaryIntegral(BottomEidM, ElementalMassMatrix2d, OP11);
        
        OP11 = LocalDownLDGBoundaryIntegral(BottomEidM, ElementalMassMatrix2d, OP11);
        %> Local Up Integral part
        OP11 = LocalUpBoundaryIntegral(UpEidM, ElementalMassMatrix2d, OP11);
        
        OP11 = LocalUpLDGBoundaryIntegral(UpEidM, ElementalMassMatrix2d, OP11);
        
        %> Assemble the local integral part into the StiffMatrix
        StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
        %> The upper adjacent cell part
        OP12 = zeros(meshUnion.cell.Np);
        OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, ElementalMassMatrix2d, OP12);
        
        OP12 = AdjacentUpLDGBoundaryIntegral(UpEidM, BottomEidM, ElementalMassMatrix2d, OP12);
        
        StiffMatrix(UpAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
        %> The lower adjacent cell part
        OP12 = zeros(meshUnion.cell.Np);
        OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, ElementalMassMatrix2d, OP12);
        
        OP12 = AdjacentDownLDGBoundaryIntegral(BottomEidM, UpEidM, ElementalMassMatrix2d, OP12);
        
        StiffMatrix(BottomAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
    end
    %> for the bottom most cell
    AdjacentRows = ((k-1)*Np*Nz + (meshUnion.Nz-2)*Np+1:(k-1)*Np*Nz + (meshUnion.Nz-1)*Np)';
    LocalRows    = ((k-1)*Np*Nz + (meshUnion.Nz-1)*Np+1:(k-1)*Np*Nz + meshUnion.Nz*Np)';
    LocalColumns = (k-1)*Np*Nz + (meshUnion.Nz-1)*Np+1:(k-1)*Np*Nz + meshUnion.Nz*Np;
    %> Volume Integral Part
    OP11 = ElementalMassMatrix3d * Dz3d;
    %> Local Up Integral part
    OP11 = LocalUpBoundaryIntegral(UpEidM, ElementalMassMatrix2d, OP11);
    
    OP11 = LocalUpLDGBoundaryIntegral(UpEidM, ElementalMassMatrix2d, OP11);
    %> Assemble the local integral part into the StiffMatrix
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP11;
    %> The upper adjacent cell part
    OP12 = zeros(meshUnion.cell.Np);
    OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, ElementalMassMatrix2d, OP12);
    
    OP12 = AdjacentUpLDGBoundaryIntegral(UpEidM, BottomEidM, ElementalMassMatrix2d, OP12);
    
    StiffMatrix(AdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix3d\OP12;
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
n0 = 6;
%> here Nz stands for ratio between area of surface and volume of the studied cell
Nz = mesh3d.Nz;
for i = 1:mesh2d.K
    %> The surface most face for each column
    Tau(1,i) = (P+1)*(P+3)/3*n0/2*(Nz + 1)*max(DiffusionCoefficient(UpEidM, (i-1)*Nz+1));
    for j = 2:Nz
        Tau(j,i) = (P+1)*(P+3)/3*n0/2*(Nz + 1)*max(max(DiffusionCoefficient(BotEidM, (i-1)*Nz+j-1)),...
            max(DiffusionCoefficient(UpEidM, (i-1)*Nz+j)));
    end
    %> The bottom most face for each column
    Tau(Nz+1,i) = (P+1)*(P+3)/3*n0/2*(Nz+1)*max(DiffusionCoefficient(BotEidM, (i-1)*Nz+Nz));
end
end

function OP11 = LocalUpBoundaryIntegral(eidM, massMatrix2d, OP11)
OP11(eidM,eidM) = OP11(eidM,eidM) -  0.5*massMatrix2d; %checked
end

function OP11 = LocalUpLDGBoundaryIntegral(eidM, massMatrix2d, OP11)
OP11(eidM,eidM) = OP11(eidM,eidM) +  0.5*massMatrix2d; %checked
end

function OP11 = LocalDownBoundaryIntegral(eidM, massMatrix2d, OP11)
OP11(eidM,eidM) = OP11(eidM,eidM) +  0.5*massMatrix2d;   %checked
end

function OP11 = LocalDownLDGBoundaryIntegral(eidM, massMatrix2d, OP11)
OP11(eidM,eidM) = OP11(eidM,eidM) +  0.5*massMatrix2d; 
end

function OP12 = AdjacentDownBoundaryIntegral(eidM, eidP, massMatrix2d, OP12)
OP12(eidP,eidM) = OP12(eidP,eidM) + 0.5*massMatrix2d;    %checked
end

function OP12 = AdjacentDownLDGBoundaryIntegral(eidM, eidP, massMatrix2d, OP12)
OP12(eidP,eidM) = OP12(eidP,eidM) - 0.5*massMatrix2d;    %checked
end

function OP12 = AdjacentUpBoundaryIntegral(eidM, eidP, massMatrix2d, OP12)
OP12(eidP,eidM) = OP12(eidP,eidM) - 0.5*massMatrix2d;          %checked
end

function OP12 = AdjacentUpLDGBoundaryIntegral(eidM, eidP, massMatrix2d, OP12)
OP12(eidP,eidM) = OP12(eidP,eidM) - 0.5*massMatrix2d;          %checked
end


function [ OP11, LBRHS ] = ImposeSurfaceDirichletBoundaryCondition(obj, UpEidM, massMatrix3d, massMatrix2d, OP11)
OP11 = DirichletTopBoundaryIntegral(UpEidM, massMatrix2d, OP11);
SurfTaux = obj.DirichExact;
temphuRHS = zeros(obj.meshUnion.cell.Np, 1);
%> move these terms about Dirichlet boundary to the right hand side, we just need to add this term to the right hand side
temphuRHS(UpEidM) = -1*massMatrix2d * SurfTaux;
LBRHS = massMatrix3d\temphuRHS;
end

function OP11 = DirichletTopBoundaryIntegral(eidM, massMatrix2d, OP11)
OP11(eidM,eidM) = OP11(eidM,eidM) - massMatrix2d; %checked
end






















