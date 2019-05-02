function fphys = matNdgConservativeNonhydrostaticUpdata(obj, physClass, fphys)
%> @brief Function to make the nonhydrostatic correction
%> @details
%> Function to make the nonhydrostatic correction
%> @param[in] physClass The hydrostatic solver
%> @param[in] fphys The fphys field
mesh = physClass.meshUnion(1);

obj.matAssembleWetDryInterface(mesh);

[UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY, UpdatedNP ]...
    = obj.matReconstructStiffmatrixRelatedMatrix( physClass);

% obj.matAssemblePointToCellInformation(mesh.K, mesh.cell.Np, UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX,...
%     UpdatedSPNPY, obj.NPBX, obj.NPBY, UpdatedFNPBX, UpdatedFNPBY, obj.NP);

StiffMatrix = obj.matAssembleConservativeGlobalSparseStiffMatrix( physClass, UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX,...
    UpdatedSPNPY, UpdatedNP, fphys);

NonhydrostaticRHS  = obj.matEvaluateConservativeNonhydrostaticRHS( fphys, physClass );

% [StiffMatrix, NonhydrostaticRHS] = obj.matResembleGlobalMatrix(mesh, StiffMatrix, NonhydrostaticRHS);

[NewStiffMatrix, NewNonhydrostaticRHS] = getStiffAndNonhydroRHS( obj, physClass, StiffMatrix, NonhydrostaticRHS );

NonhydrostaticPressure = obj.matCalculateNonhydrostaticPressure(NewStiffMatrix, NewNonhydrostaticRHS);
fphys = obj.matUpdateConservativeFinalVelocity( NonhydrostaticPressure , physClass, fphys);
obj.ZeroFluxBoundary = []; obj.ZeroFluxBoundaryIndex = 0;
obj.AdjacentDryCellAndFace = [];obj.NonhydroFmPoint = [];
obj.NonhydroFpPoint = [];
end

function [NewStiffMatrix, NewNonhydrostaticRHS] = getStiffAndNonhydroRHS( obj, physClass, StiffMatrix, NonhydrostaticRHS )
% VolumeIntegral = zeros(size(obj.Vq{1},2), size(Variable,2));
Np = physClass.meshUnion(1).cell.Np;
NewStiffMatrix = zeros(numel(physClass.meshUnion(1).x));
for i = 1:size(StiffMatrix,2)
%     display(i);
    tempStiffMatrix = reshape(StiffMatrix(:,i), size(obj.Vq{1}, 1), physClass.meshUnion(1).K);
    for j = 1:physClass.meshUnion(1).K
        data = ( obj.Vq{1} )' * diag(obj.wJ{1}(:,j)) * tempStiffMatrix(:,j) ;
        NewStiffMatrix((j-1)*Np+1 : j *Np,i) = data(:);
    end
end

% tempNonhydrostaticRHS = reshape(NonhydrostaticRHS, size(obj.Vq{1}, 1), physClass.meshUnion(1).K);
NewNonhydrostaticRHS =  zeros(numel(physClass.meshUnion(1).x), 1);
for k = 1:physClass.meshUnion(1).K
    data = ( obj.Vq{1} )' * diag(obj.wJ{1}(:,k)) * NonhydrostaticRHS(:,k);
    NewNonhydrostaticRHS((k-1)*Np+1 : k *Np)  = data(:);
end

end