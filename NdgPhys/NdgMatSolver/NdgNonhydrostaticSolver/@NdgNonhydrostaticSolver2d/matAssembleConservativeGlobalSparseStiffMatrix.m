function StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, fphys)
%> @brief Function to assemble the global stiff matrix
%> @details
%> Function to assemble the global sparse stiff matrix
%> @param[in] UpdatedPNPX The first order derivative of Nonhydrostatic pressure with regard to x
%> @param[in] UpdatedPNPY The first order derivative of Nonhydrostatic pressure with regard to y
%> @param[in] UpdatedSPNPX The second order derivative of Nonhydrostatic pressure with regard to x
%> @param[in] UpdatedSPNPY The second order derivative of Nonhydrostatic pressure with regard to y
%> @param[in] UpdatedFNPBX The product of Nonhydrostatic pressure with bottom gradient in x direction and this term is included in a flux term
%> @param[in] UpdatedFNPBY The product of Nonhydrostatic pressure with bottom gradient in y direction and this term is included in a flux term
%> @param[in] fphys The physical field
%> @param[in] PhysClass The hydrostatic Solver
% mesh = PhysClass.meshUnion(1);

h = fphys{1}(:,:,1);

% StiffMatrix = mxAssembleGlobalStiffMatrix(obj.dt, h, UpdatedSPNPX, UpdatedSPNPY, UpdatedPNPX, obj.fhx, UpdatedPNPY, obj.fhy,...
%     obj.NP, obj.H2Bx, obj.H2By, ( obj.HBx ).^2, ( obj.HBy ).^2,  obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix);

StiffMatrix = mxAssembleGlobalStiffMatrix(obj.dt, h, obj.TempSecondOrderTerm, obj.TempPNPX, obj.fhx, obj.TempPNPY, obj.fhy,...
    obj.NP, obj.H2Bx, obj.H2By, ( obj.HBx ).^2, ( obj.HBy ).^2,  obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix);

end