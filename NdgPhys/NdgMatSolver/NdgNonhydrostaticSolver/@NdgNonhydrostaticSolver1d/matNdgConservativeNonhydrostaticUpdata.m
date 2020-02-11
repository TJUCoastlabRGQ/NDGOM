function fphys = matNdgConservativeNonhydrostaticUpdata(obj, physClass, fphys) 
%> @brief Function to make the nonhydrostatic correction
%> @details
%> Function to make the nonhydrostatic correction
%> @param[in] physClass The hydrostatic solver
%> @param[in] fphys The fphys field
mesh = physClass.meshUnion(1);

obj.matAssembleWetDryInterface1d(mesh);    

obj.matReconstructStiffmatrixRelatedMatrix( mesh ); 

obj.matAssemblePointToCellInformation(mesh.K, mesh.cell.Np); 

obj.matCalculateFphysDerivative( mesh, fphys, physClass); 

StiffMatrix = obj.matAssembleConservativeGlobalSparseStiffMatrix( fphys ); 

NonhydrostaticRHS  = obj.matEvaluateConservativeNonhydrostaticRHS( fphys, physClass ); 

[StiffMatrix, NonhydrostaticRHS] = obj.matResembleGlobalMatrix(mesh, StiffMatrix, NonhydrostaticRHS);

NonhydrostaticPressure = obj.matCalculateNonhydrostaticPressure(StiffMatrix, NonhydrostaticRHS);

fphys = obj.matUpdateConservativeFinalVelocity( NonhydrostaticPressure , physClass, fphys); 


end