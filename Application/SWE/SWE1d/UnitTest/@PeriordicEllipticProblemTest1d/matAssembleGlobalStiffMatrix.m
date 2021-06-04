function [ StiffMatrix, LRHS, RRHS ] = matAssembleGlobalStiffMatrix(obj)
meshUnion = obj.meshUnion;
LeftEidM = meshUnion.cell.Fmask(1);
RightEidM= meshUnion.cell.Fmask(2);
%> Calculation of the penalty parameter
Tau = zeros(1, meshUnion.K);
DiffusionCoefficient  = ones(size(meshUnion.x));
Tau = matCalculatePenaltyParameter(meshUnion, DiffusionCoefficient, LeftEidM, RightEidM, Tau);
Nz = meshUnion.K;
Np = meshUnion.cell.Np;
StiffMatrix     = zeros(Np*Nz);
%> At present, we assume the mesh is uniform in the vertical direction
ElementalMassMatrix = diag(meshUnion.J(:,1))*meshUnion.cell.M;
ElementalMassMatrix2d = 1;
Dz1d = diag(meshUnion.rx(:,1))*meshUnion.cell.Dr;
%> first cell first, then the other cells left
LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,1)) * Dz1d;
AdjacentPhysicalDiffMatrix = diag(DiffusionCoefficient(:,2)) * Dz1d;
UpAdjacentRows = ((meshUnion.K-1)*Np+1:meshUnion.K*Np)';
LocalRows    = (1:Np)';
LocalColumns = 1:Np;
BottomAdjacentRows = (Np+1:2*Np)';
%> Volume Integral Part, $$-\int_{\Omega}\nabla_hs\cdot\nabla_h p_hdx$$
OP11 = -Dz1d' * ElementalMassMatrix * LocalPhysicalDiffMatrix;
%> Local Bottom Integral part, checked
% [LRHS] = ImposeSurfaceNewmannBoundaryCondition( obj, LeftEidM, ElementalMassMatrix2d, ElementalMassMatrix );

[ OP11, LRHS ] = ImposeSurfaceDirichletBoundaryCondition(obj, LeftEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, ElementalMassMatrix, Tau(1), OP11);

%> Local Bottom Integral part
OP11 = LocalRightBoundaryIntegral(RightEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(2), OP11);

StiffMatrix(LocalRows(:),LocalColumns(:)) =  ElementalMassMatrix\OP11;
%> Adjacent bottom integral part
OP12 = zeros(meshUnion.cell.Np);
OP12 = AdjacentRightBoundaryIntegral(RightEidM, LeftEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP12);
StiffMatrix(BottomAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
% OP12 = zeros(meshUnion.cell.Np);
% OP12 = AdjacentLeftBoundaryIntegral(LeftEidM, RightEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP12);
% StiffMatrix(UpAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
for j = 2:meshUnion.K-1
    UpAdjacentRows = ((j-2)*Np+1:(j-1)*Np)';
    LocalRows    = ((j-1)*Np+1:j*Np)';
    BottomAdjacentRows = (j*Np+1:(j+1)*Np)';
    LocalColumns = (j-1)*Np+1:j*Np;
    UpPhysicalDiffMatrix = diag(DiffusionCoefficient(:,j-1)) * Dz1d;
    LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,j)) * Dz1d;
    BottomPhysicalDiffMatrix = diag(DiffusionCoefficient(:,j+1)) * Dz1d;
    %> Volume Integral Part, $$-\int_{\Omega}\nabla_hs\cdot\nabla_h p_hdx$$ 
    OP11 = -Dz1d' * ElementalMassMatrix * LocalPhysicalDiffMatrix;
    %> Local Bottom Integral part, checked
    OP11 = LocalLeftBoundaryIntegral(LeftEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP11);
    %> Local Up Integral part, checked
    OP11 = LocalRightBoundaryIntegral(RightEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP11);
    %> Assemble the local integral part into the StiffMatrix
    StiffMatrix(LocalRows(:),LocalColumns(:)) = ElementalMassMatrix\OP11;
    %> The upper adjacent cell part
    OP12 = zeros(meshUnion.cell.Np);
    %> Influence of local element have on right adjacent cell , checked.
    OP12 = AdjacentRightBoundaryIntegral(RightEidM, LeftEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP12);
    StiffMatrix(BottomAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
    %> The lower adjacent cell part
    OP12 = zeros(meshUnion.cell.Np);
    %> Influence of local element have on left adjacent cell, checked
    OP12 = AdjacentLeftBoundaryIntegral(LeftEidM, RightEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP12);
    StiffMatrix(UpAdjacentRows(:),LocalColumns(:)) =  ElementalMassMatrix\OP12;
end
%> for the bottommost(x = 0) cell
UpAdjacentRows = ((meshUnion.K-2)*Np+1:(meshUnion.K-1)*Np)';
LocalRows    = ((meshUnion.K-1)*Np + 1:meshUnion.K * Np)';
LocalColumns = (meshUnion.K-1)*Np + 1:meshUnion.K * Np;

%> Volume Integral Part, $$-\int_{\Omega}\nabla_hs\cdot\nabla_h p_hdx$$
OP11 = -Dz1d' * ElementalMassMatrix * LocalPhysicalDiffMatrix;
%> Local Bottom Integral part, checked
OP11 = LocalLeftBoundaryIntegral(LeftEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(1), OP11);

% [RRHS] = ImposeBottomNewmannBoundaryCondition( obj, RightEidM, ElementalMassMatrix2d, ElementalMassMatrix );

[ OP11, RRHS ] = ImposeBottomDirichletBoundaryCondition(obj, RightEidM, LocalPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, ElementalMassMatrix, Tau(end), OP11);

%> Local Bottom Integral part
StiffMatrix(LocalRows(:),LocalColumns(:)) =  ElementalMassMatrix\OP11;
%> Adjacent bottom integral part
% OP12 = zeros(meshUnion.cell.Np);
% OP12 = AdjacentRightBoundaryIntegral(RightEidM, LeftEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(2), OP12);
% StiffMatrix(BottomAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
OP12 = zeros(meshUnion.cell.Np);
OP12 = AdjacentLeftBoundaryIntegral(LeftEidM, RightEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz1d, ElementalMassMatrix2d, Tau(2), OP12);
StiffMatrix(UpAdjacentRows(:),LocalColumns(:)) = ElementalMassMatrix\OP12;
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

for i = 1:numel(Tau)
    Tau(i) = (P+1)*(P+1)/1*n0/2*1/mesh.LAV(1)*DiffusionCoefficient(1);
end

end

function OP11 = LocalRightBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
% $$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$$
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * 1 * 0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
% $$\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$$
OP11(eidM, :)   = OP11(eidM, :)   +  0.5*massMatrix2d*physicalDiffMatrix(eidM,:); %checked
% $$\int_{\epsilon_i}\tau^k[p_h][s]d\boldsymbol{x}$$
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d; %checked
end

function OP11 = LocalLeftBoundaryIntegral(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
% $$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$$
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * (-1) *  0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
% $$\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$$
OP11(eidM, :)   = OP11(eidM, :)   -  0.5*massMatrix2d*physicalDiffMatrix(eidM,:);  %checked
% $$\int_{\epsilon_i}\tau^k[p_h][s]d\boldsymbol{x}$$
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d;   %checked
end

function OP12 = AdjacentLeftBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
%> Here, Down or up is relative to local cell
epsilon = -1;
% $$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$$
OP12(:,eidM)    = OP12(:,eidM) - epsilon * (-1) * 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;
% $$\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$$
OP12(eidP,:)    = OP12(eidP,:) +  0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);  %checked
% $$\int_{\epsilon_i}\tau^k[p_h][s]d\boldsymbol{x}$$
OP12(eidP,eidM) = OP12(eidP,eidM) +  Tau * massMatrix2d;    %checked
end

function OP12 = AdjacentRightBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, Dz, massMatrix2d, Tau, OP12)
epsilon = -1;
%$$\int_{\epsilon_i}[p_h]\cdot\left\{\nabla_h s\right\}d\boldsymbol{x}$$
OP12(:,eidM)    = OP12(:,eidM) -  epsilon * (1) * 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;   %checked
% $$\int_{\epsilon_i}\left\{\nabla_hp_h\right\}[s]d\boldsymbol{x}$$
OP12(eidP,:)    = OP12(eidP,:) -  0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);    %checked
% $$\int_{\epsilon_i}\tau^k[p_h][s]d\boldsymbol{x}$$
OP12(eidP,eidM) = OP12(eidP,eidM) +  Tau * massMatrix2d;          %checked
end

function [ OP11, RBRHS ] = ImposeBottomDirichletBoundaryCondition(obj, RightEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, massMatrix, Tau, OP11)
epsilon = -1;
OP11 = DirichletRightBoundaryCondition(RightEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
BottomValue = obj.DirichletExact;
TempRBRHS = zeros(obj.meshUnion.cell.Np, 1);
TempRBRHS(RightEidM) = -1*Tau*massMatrix2d * BottomValue;
TempRBRHS = TempRBRHS + epsilon * LocalPhysicalDiffMatrix(RightEidM,:)'*massMatrix2d*(-1*BottomValue);
RBRHS = massMatrix\TempRBRHS;
end

function [ OP11, LBRHS ] = ImposeSurfaceDirichletBoundaryCondition(obj, LeftEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, massMatrix, Tau, OP11)
epsilon = -1;
OP11 = DirichletLeftBoundaryCondition(LeftEidM, LocalPhysicalDiffMatrix, Dz3d, massMatrix2d, Tau, OP11);
SurfaceValue = obj.DirichletLeftExact;
TempRBRHS = zeros(obj.meshUnion.cell.Np, 1);
%> move these terms about Dirichlet boundary to the left 
TempRBRHS(LeftEidM) = 1*Tau*massMatrix2d * SurfaceValue;
TempRBRHS = TempRBRHS + epsilon * LocalPhysicalDiffMatrix(LeftEidM,:)'*massMatrix2d*(-1*SurfaceValue);
LBRHS = massMatrix\TempRBRHS;
end

function OP11 = DirichletLeftBoundaryCondition(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
%$$\int_{\partial \Omega^D}p_h\nabla_hs\cdot \boldsymbol{n}d\boldsymbol{x}$$
OP11(:, eidM)   = OP11(:, eidM)   + epsilon * 1 * physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
%$$\int_{\partial\Omega^D}s\boldsymbol{n}\cdot\nabla_hp_hd\boldsymbol{x}$$
OP11(eidM, :)   = OP11(eidM, :)   -  massMatrix2d*physicalDiffMatrix(eidM,:); %checked
%$$\int_{\partial \Omega^D}\tau^ksp_hd\boldsymbol{x}$$
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d; %checked
end

function [LBRHS] = ImposeSurfaceNewmannBoundaryCondition( obj, eidM, massMatrix2d, massMatrix )
%> This part is negative, as this is teated explicitly
WindTaux = obj.NewmannExact;
temphuRHS = zeros(obj.meshUnion.cell.Np, 1);
temphuRHS(eidM) = massMatrix2d * WindTaux;
LBRHS = massMatrix\temphuRHS;
% LBRHS = temphuRHS;
end

function [RBRHS] = ImposeBottomNewmannBoundaryCondition( obj, eidM, massMatrix2d, massMatrix )
%> This part is negative, as this is teated explicitly
WindTaux = obj.NewmannRightExact;
temphuRHS = zeros(obj.meshUnion.cell.Np, 1);
temphuRHS(eidM) = massMatrix2d * WindTaux;
RBRHS = massMatrix\temphuRHS;
RBRHS = -1 * RBRHS;
end

function OP11 = DirichletRightBoundaryCondition(eidM, physicalDiffMatrix, Dz, massMatrix2d, Tau, OP11)
epsilon = -1;
%$$\int_{\partial \Omega^D}p_h\nabla_hs\cdot \boldsymbol{n}d\boldsymbol{x}$$
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * 1 * physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
%$$\int_{\partial\Omega^D}s\boldsymbol{n}\cdot\nabla_hp_hd\boldsymbol{x}$$
OP11(eidM, :)   = OP11(eidM, :)   +  massMatrix2d*physicalDiffMatrix(eidM,:); %checked
%$$\int_{\partial \Omega^D}\tau^ksp_hd\boldsymbol{x}$$
OP11(eidM,eidM) = OP11(eidM,eidM) -  Tau*massMatrix2d; %checked
end

