function fphys = matCalculateImplicitRHS( obj, physClass, DiffusionCoefficient, SystemRHS, ImplicitParameter, dt, RKIndex, IMStage)
%> @brief Calculating the right hand side corresponding to the vertical diffusion term and
%> return the physical field with vertical diffusion considered
%> @detail this function is used to calculate the right hand side corresponding to the vertical
%> diffusion term and return the updated physical field at each Runge-Kutta time stage
%> @param[in] physClass The physical solver establised
%> @param[in] Height The water depth
%> @param[in] SystemRHS The right hand side corresponding to the discretization of vertical diffusion term
%> @param[in] ImplicitParameter The implicit parameter at the corresponding IMEXRK stage
%> @param[in] dt The time step
%> @param[in] RKIndex The intermediate stage of IMEXRK time stepping method
%> @param[in] IMStage The total stage of the implicit part of IMEXRK time stepping method
%> @param[out] fphys The physical field with vertical diffusion
%> considered

BottomEidM   = physClass.meshUnion(1).cell.Fmask(physClass.meshUnion(1).cell.Fmask(:,end-1)~=0,end-1);
UpEidM     = physClass.meshUnion(1).cell.Fmask(physClass.meshUnion(1).cell.Fmask(:,end)~=0,end);
K = physClass.meshUnion(1).K;
Nz = physClass.meshUnion(1).Nz;
Np = physClass.meshUnion(1).cell.Np;
% ImplicithuRHS3d = zeros(size(meshUnion.x));
% ImplicithvRHS3d = zeros(size(meshUnion.x));
% hu              = zeros(size(meshUnion.x));
% hv              = zeros(size(meshUnion.x));
fphys = zeros( Np, physClass.meshUnion(1).K, physClass.Nvar );
StiffMatrix     = zeros( Np*Nz, Np*Nz, physClass.Nvar );
for i =1:physClass.meshUnion(1).mesh2d(1).K
    %> At present, we assume the mesh is uniform in the vertical direction
    ElementalMassMatrix3d = diag(physClass.meshUnion(1).J(:,(i-1)*Nz+1)) * physClass.meshUnion(1).cell.M;
    ElementalMassMatrix2d = diag(physClass.meshUnion(1).mesh2d(1).J(:,i)) * physClass.meshUnion(1).mesh2d(1).cell.M;
    Dz3d = diag(physClass.meshUnion(1).tz(:,(i-1)*Nz+1)) * physClass.meshUnion(1).cell.Dt;
    %> first cell first, then the other cells left
    LocalRows    = (1:Np)';
    AdjacentRows = (Np+1:2*Np)';
    LocalColumns = 1:Np;
    LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+1))*Dz3d;
    AdjacentPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+2))*Dz3d;
    %> Volume Integral Part
    OP11 = -Dz3d' * ElementalMassMatrix3d * LocalPhysicalDiffMatrix;
    %> Local Bottom Integral part
%     obj.tau = zeros( numel(BotEidM), physClass.mesh2d(1).K * ( physClass.meshUnion(1).Nz+1 ) )
%     OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(2,i), OP11);
    OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + 2), OP11);
    %> Adjacent bottom integral part
    OP12 = zeros(physClass.meshUnion(1).cell.Np);
%     OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(2,i), OP12);
    OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + 2), OP12);
    %> Boundary part, not considered here
    [ SystemRHS(:,(i-1) * Nz + 1,:), SurfStiffMatrix ] = ImposeNewmannBoundaryCondition(UpEidM, physClass.SurfBoundNewmannDate(:,i,:),...
        ElementalMassMatrix2d, ElementalMassMatrix3d, dt, ImplicitParameter, SystemRHS(:,(i-1) * Nz + 1,:));
    %> Assemble into the StiffMatrix, it's noted that, the diagonal part has been included by eye(Np)
    for var = 1:2
        StiffMatrix(LocalRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\OP11;
        StiffMatrix(AdjacentRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\OP12;
    end
    for var = 3:physClass.Nvar
        StiffMatrix(LocalRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\(OP11./obj.Prantl);
        StiffMatrix(AdjacentRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\(OP12./obj.Prantl);
    end
    for j = 2:physClass.meshUnion(1).Nz-1
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
%         OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(j+1,i), OP11);
        OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + j + 1 ), OP11);
        %> Local Up Integral part
%         OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(j,i), OP11);
        OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + j ), OP11);

        %> Assemble the local integral part into the StiffMatrix
        %         StiffMatrix(LocalRows(:),LocalColumns(:),1:2) = ElementalMassMatrix3d\OP11;
        %         StiffMatrix(LocalRows(:),LocalColumns(:),3:end) = ElementalMassMatrix3d\(OP11./obj.Prantl);
        for var = 1:2
            StiffMatrix(LocalRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\OP11;
        end
        for var = 3:physClass.Nvar
            StiffMatrix(LocalRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\(OP11./obj.Prantl);
        end
        
        %> The upper adjacent cell part
        OP12 = zeros(physClass.meshUnion(1).cell.Np);
%         OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(j,i), OP12);
        OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, LocalPhysicalDiffMatrix, UpPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + j ), OP12);
        
        for var = 1:2
            StiffMatrix(UpAdjacentRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\OP12;
        end
        for var = 3:physClass.Nvar
            StiffMatrix(UpAdjacentRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\(OP12./obj.Prantl);
        end
        
        %> The lower adjacent cell part
        OP12 = zeros(Np);
%         OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(j+1,i), OP12);
        OP12 = AdjacentDownBoundaryIntegral(BottomEidM, UpEidM, LocalPhysicalDiffMatrix, BottomPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + j + 1 ), OP12);

        for var = 1:2
            StiffMatrix(BottomAdjacentRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\OP12;
        end
        for var = 3:physClass.Nvar
            StiffMatrix(BottomAdjacentRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\(OP12./obj.Prantl);
        end
    end
    %> for the bottom most cell
    AdjacentRows = ((Nz-2)*Np+1:(Nz-1)*Np)';
    LocalRows    = ((Nz-1)*Np+1:Nz*Np)';
    LocalColumns = (Nz-1)*Np+1:Nz*Np;
    LocalPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+Nz))*Dz3d;
    AdjacentPhysicalDiffMatrix = diag(DiffusionCoefficient(:,(i-1)*Nz+Nz-1))*Dz3d;
    %> Volume Integral Part
    OP11 = -Dz3d' * ElementalMassMatrix3d * LocalPhysicalDiffMatrix;
    %> Local Up Integral part
%     OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(Nz,i), OP11);
    OP11 = LocalUpBoundaryIntegral(UpEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + Nz ), OP11);
    
    %> Impose Newmann boundary for the third to the last physical field
    for var = 3:physClass.Nvar
        StiffMatrix(LocalRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\(OP11./obj.Prantl);
    end
    [ SystemRHS(:,i*Nz,:), BotStiffMatrix ] = ImposeNewmannBoundaryCondition(BottomEidM, physClass.BotBoundNewmannDate(:,i,:),...
        ElementalMassMatrix2d, ElementalMassMatrix3d, dt, ImplicitParameter, SystemRHS(:,i*Nz,:));
    %% This part is used to impose Newmann boundary condition for hu and hv
%     for var = 1:2
%         StiffMatrix(LocalRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\OP11;
%     end
    %% This part is used to impose homogeneous Dirichlet boundary condition for hu and hv
    %> Impose bottom boundary condition
%     [ OP11 ] = ImposeBottomDirichletBoundaryCondition(BottomEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(Nz+1,i), OP11);
    [ OP11 ] = ImposeBottomDirichletBoundaryCondition(BottomEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + Nz + 1), OP11);

    %> Assemble the local integral part into the StiffMatrix
    for var = 1:2
        StiffMatrix(LocalRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\OP11;
    end
    %%
    
    %> The upper adjacent cell part
    OP12 = zeros(Np);
%     OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(Nz,i), OP12);
    OP12 = AdjacentUpBoundaryIntegral(UpEidM, BottomEidM, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, ElementalMassMatrix2d, obj.tau(:,(i-1)*( physClass.meshUnion(1).Nz+1 ) + Nz), OP12);
    
    for var = 1:2
        StiffMatrix(AdjacentRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\OP12;
    end
    for var = 3:physClass.Nvar
        StiffMatrix(AdjacentRows(:),LocalColumns(:),var) = ElementalMassMatrix3d\(OP12./obj.Prantl);
    end
    %
    for var = 1:physClass.Nvar
        temphuRHS = SystemRHS(:,(i-1)*Nz + 1 : i*Nz,var);
        fphys( (var-1)* K * Np + (i-1)*Nz * Np + 1  : (var-1)* K * Np + i*Nz * Np)...
            = (eye(Np*Nz) - ImplicitParameter * dt * StiffMatrix(:,:,var))\temphuRHS(:);
        physClass.ImplicitRHS( ((var-1)* ( IMStage - 1 ) + RKIndex - 1 )* K * Np + (i-1)*Nz * Np + 1 :...
            ( (var-1) * ( IMStage - 1) + RKIndex - 1 )* K * Np + i*Nz * Np ) = ...
            StiffMatrix(:,:,var) * fphys( (var-1)* K * Np + (i-1)*Nz * Np + 1  : (var-1)* K * Np + i*Nz * Np)';
        
        physClass.ImplicitRHS(:,(i-1)*Nz + 1,(var-1) * ( IMStage - 1 ) + RKIndex ) = physClass.ImplicitRHS(:,(i-1)*Nz + 1,(var-1)*( IMStage - 1 ) + RKIndex ) + SurfStiffMatrix(:,:,var);
        physClass.ImplicitRHS(:,i*Nz,(var-1)* ( IMStage - 1 )  + RKIndex ) = physClass.ImplicitRHS(:,i*Nz,(var-1)*( IMStage - 1 ) + RKIndex ) + BotStiffMatrix(:,:,var);
    end
    
end
end

function OP11 = LocalUpBoundaryIntegral(eidM, physicalDiffMatrix, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * 1 * 0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   + 0.5*massMatrix2d*physicalDiffMatrix(eidM,:); %checked
OP11(eidM,eidM) = OP11(eidM,eidM) - diag(Tau)*massMatrix2d; %checked
end

function OP11 = LocalDownBoundaryIntegral(eidM, physicalDiffMatrix, massMatrix2d, Tau, OP11)
epsilon = -1;
OP11(:, eidM)   = OP11(:, eidM)   - epsilon * (-1) *  0.5*physicalDiffMatrix(eidM,:)'*massMatrix2d; %checked
OP11(eidM, :)   = OP11(eidM, :)   -  0.5*massMatrix2d*physicalDiffMatrix(eidM,:);  %checked
OP11(eidM,eidM) = OP11(eidM,eidM) -  diag(Tau)*massMatrix2d;   %checked
end

function OP12 = AdjacentDownBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, massMatrix2d, Tau, OP12)
%> Here, Down or up is relative to local cell
epsilon = -1;
OP12(:,eidM)    = OP12(:,eidM) - epsilon * (-1) * 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;
OP12(eidP,:)    = OP12(eidP,:) +  0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);  %checked
OP12(eidP,eidM) = OP12(eidP,eidM) +  diag(Tau) * massMatrix2d;    %checked
end

function OP12 = AdjacentUpBoundaryIntegral(eidM, eidP, LocalPhysicalDiffMatrix, AdjacentPhysicalDiffMatrix, massMatrix2d, Tau, OP12)
epsilon = -1;
OP12(:,eidM)    = OP12(:,eidM) -  epsilon * (1) * 0.5 * AdjacentPhysicalDiffMatrix(eidP,:)'*massMatrix2d;   %checked
OP12(eidP,:)    = OP12(eidP,:) -  0.5 * massMatrix2d * LocalPhysicalDiffMatrix(eidM,:);    %checked
OP12(eidP,eidM) = OP12(eidP,eidM) +  diag(Tau) * massMatrix2d;          %checked
end


function [ huRHS, huStiffMatrix ] = ImposeNewmannBoundaryCondition(eidM, BoundNewmannDate, massMatrix2d, massMatrix3d, dt, ImplicitParameter, huRHS)
%> This part is negative, as this is teated explicitly
huStiffMatrix = zeros( size( huRHS ) );
for i = 1:size(huRHS,3)
    temphuRHS = zeros(size(huRHS(:,:,i)));
    temphuRHS(eidM) = massMatrix2d * BoundNewmannDate(:,:,i);
    huStiffMatrix(:,:,i) = massMatrix3d\temphuRHS;
    huRHS(:,:,i) = huRHS(:,:,i) + dt*ImplicitParameter*huStiffMatrix(:,:,i);
end
end

function [ OP11 ] = ImposeBottomDirichletBoundaryCondition(BottomEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, Tau, OP11)
OP11 = LocalDownBoundaryIntegral(BottomEidM, LocalPhysicalDiffMatrix, ElementalMassMatrix2d, 2 * Tau, OP11);
end