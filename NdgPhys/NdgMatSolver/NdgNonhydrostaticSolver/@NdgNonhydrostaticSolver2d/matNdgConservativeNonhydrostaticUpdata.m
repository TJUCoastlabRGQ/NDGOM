function fphys = matNdgConservativeNonhydrostaticUpdata(obj, physClass, fphys)   %protected function in nonhydrostait
%> @brief Function to make the nonhydrostatic correction
%> @details
%> Function to make the nonhydrostatic correction
%> @param[in] physClass The hydrostatic solver
%> @param[in] fphys The fphys field
mesh = physClass.meshUnion(1);

obj.matAssembleWetDryInterface2d(mesh);    % move to 1d and 2d solver, actually this part is same for 1d and 2d solver

obj.matReconstructStiffmatrixRelatedMatrix;  % abstract in abstract nonhydrostatic solver,  different for Gauss Quad and Quad Free

obj.matAssemblePointToCellInformation(mesh.K, mesh.cell.Np);  % protected for abstract nonhydrostatic solver, reload in quad free version

obj.matCalculateFphysDerivative( mesh, fphys, physClass);  % abstract in  abstract nonhydrostatic solver, 

StiffMatrix = obj.matAssembleConservativeGlobalSparseStiffMatrix( fphys ); % abstract in abstract nonhydrostatic solver, different for Gauss quad and quad free version

NonhydrostaticRHS  = obj.matEvaluateConservativeNonhydrostaticRHS( fphys, physClass ); % abstract in abstract nonhydrostatic solver, different for Gauss quad and quad free version

[StiffMatrix, NonhydrostaticRHS] = obj.matResembleGlobalMatrix(mesh, StiffMatrix, NonhydrostaticRHS);% move to abstract nonhydrostatic solver, the same for 1d and 2d solver

NonhydrostaticPressure = obj.matCalculateNonhydrostaticPressure(StiffMatrix, NonhydrostaticRHS);% move to abstract nonhydrostatic solver, the same for 1d and 2d solver
fphys = obj.matUpdateConservativeFinalVelocity( NonhydrostaticPressure , physClass, fphys); % abstract  in abstract nonhydrostatic solver, different for Gauss quad and quad free version


end