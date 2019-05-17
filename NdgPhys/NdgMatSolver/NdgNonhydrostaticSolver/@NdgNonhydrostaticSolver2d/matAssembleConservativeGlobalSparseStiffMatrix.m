function StiffMatrix = matAssembleConservativeGlobalSparseStiffMatrix(obj, UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY, fphys, PhysClass)
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
mesh = PhysClass.meshUnion(1);

% rho = PhysClass.rho;
h = fphys{1}(:,:,1);
BoundaryEdge = PhysClass.meshUnion(1).BoundaryEdge;
InnerEdge = PhysClass.meshUnion(1).InnerEdge;

[ fhx, fhy ] = obj.matCalculateConservativeVariableRHSMatrix( PhysClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 1);

[ obj.HBx, obj.HBy ] = obj.matCalculateCharacteristicMatrix( mesh,  BoundaryEdge, InnerEdge, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]),  num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), enumNonhydroBoundaryCondition.Zero);  

% [TempHBx, TempHBy] = obj.matCalculateLDGAuxialaryVariable( mesh, BoundaryEdge, InnerEdge, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]));
% [ H2Bx, H2By ] = obj.matCalculateLDGSecondOrderVariable( mesh, BoundaryEdge, InnerEdge, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), num2cell(TempHBx,[1 2]), num2cell(TempHBy,[1 2]) );

[ TempHBx, TempHBy ]  = obj.matCalculateCharacteristicMatrix( mesh, BoundaryEdge, InnerEdge, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), enumNonhydroBoundaryCondition.Zero);
[ H2Bx, H2By ]  = obj.matCalculateCharacteristicMatrix( mesh, BoundaryEdge, InnerEdge, num2cell(TempHBx,[1 2]), num2cell(TempHBy,[1 2]), enumNonhydroBoundaryCondition.ZeroGrad);



% obj.NP =  obj.NP;
% StiffMatrix = mxAssembleGlobalStiffMatrix(obj.dt, rho, obj.NP, h, UpdatedPNPX, obj.bx,...
%     UpdatedPNPY, obj.by, UpdatedSPNPX, UpdatedSPNPY, UpdatedFNPBX, UpdatedFNPBY, ...
%     fhx, fhy, obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix);

StiffMatrix = mxAssembleGlobalStiffMatrix(obj.dt, h, UpdatedSPNPX, UpdatedSPNPY, UpdatedPNPX, fhx, UpdatedPNPY, fhy,...
    obj.NP, H2Bx, H2By, ( obj.HBx ).^2, ( obj.HBy ).^2,  obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix);
end