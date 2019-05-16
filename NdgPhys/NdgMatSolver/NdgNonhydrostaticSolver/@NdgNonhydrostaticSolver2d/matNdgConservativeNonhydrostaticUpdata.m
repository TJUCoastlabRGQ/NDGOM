function fphys = matNdgConservativeNonhydrostaticUpdata(obj, physClass, fphys)
%> @brief Function to make the nonhydrostatic correction
%> @details
%> Function to make the nonhydrostatic correction
%> @param[in] physClass The hydrostatic solver
%> @param[in] fphys The fphys field
mesh = physClass.meshUnion(1);

obj.matAssembleWetDryInterface(mesh);

[UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY]...
    = obj.matReconstructStiffmatrixRelatedMatrix( physClass);

obj.matAssemblePointToCellInformation(mesh.K, mesh.cell.Np, UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX,...
    UpdatedSPNPY, obj.NPBX, obj.NPBY, UpdatedFNPBX, UpdatedFNPBY, obj.NP);

StiffMatrix = obj.matAssembleConservativeGlobalSparseStiffMatrix(UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX,...
    UpdatedSPNPY, UpdatedFNPBX, UpdatedFNPBY, fphys, physClass);

NonhydrostaticRHS  = obj.matEvaluateConservativeNonhydrostaticRHS( fphys, physClass );

[StiffMatrix, NonhydrostaticRHS] = obj.matResembleGlobalMatrix(mesh, StiffMatrix, NonhydrostaticRHS);

NonhydrostaticPressure = obj.matCalculateNonhydrostaticPressure(StiffMatrix, NonhydrostaticRHS);
fphys = obj.matUpdateConservativeFinalVelocity( NonhydrostaticPressure , physClass, fphys);
obj.ZeroFluxBoundary = []; obj.ZeroFluxBoundaryIndex = 0;
obj.AdjacentDryCellAndFace = [];obj.NonhydroFmPoint = [];
obj.NonhydroFpPoint = [];
end