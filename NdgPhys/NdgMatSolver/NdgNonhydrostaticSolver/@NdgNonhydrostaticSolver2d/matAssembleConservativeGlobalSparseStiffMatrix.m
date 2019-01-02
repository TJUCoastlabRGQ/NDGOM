function StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY, UpdatedFNPBX, UpdatedFNPBY...
    , fphys, PhysClass)
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
rho = PhysClass.rho;
h = fphys{1}(:,:,1);
% BoundaryEdge = PhysClass.meshUnion(1).BoundaryEdge;
% InnerEdge = PhysClass.meshUnion(1).InnerEdge;
% fhx = obj.matCalculateConservativeVariableRelatedMatrixX( PhysClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 1);
% fhy = obj.matCalculateConservativeVariableRelatedMatrixY( PhysClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 1);
[fhx, fhy] = obj.matEvaluateLocalDerivativeTerm( PhysClass.meshUnion(1), h );
StiffMatrix = mxAssembleGlobalStiffMatrix(obj.dt, rho, obj.NP, h, UpdatedPNPX, obj.bx,...
    UpdatedPNPY, obj.by, UpdatedSPNPX, UpdatedSPNPY, UpdatedFNPBX, UpdatedFNPBY, ...
    fhx, fhy, obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix);
end