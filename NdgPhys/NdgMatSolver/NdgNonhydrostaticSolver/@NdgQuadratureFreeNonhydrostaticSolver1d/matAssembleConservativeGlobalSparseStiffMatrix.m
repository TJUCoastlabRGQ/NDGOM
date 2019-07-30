function StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, fphys)
%> @brief Function to assemble the global stiff matrix
%> @details
%> Function to assemble the global sparse stiff matrix
%> @param[in] fphys The physical field
%> @param[in] PhysClass The hydrostatic Solver
%> @param[in] StiffMatrix The global stiff matrix
% mesh = PhysClass.meshUnion(1);

h = fphys{1}(:,:,1);

StiffMatrix = mxAssembleGlobalStiffMatrix(obj.dt, h, obj.TempSecondOrderTerm, obj.TempPNPX, obj.fhx, ...
    obj.NP, obj.H2Bx, ( obj.HBx ).^2,  obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix);

end