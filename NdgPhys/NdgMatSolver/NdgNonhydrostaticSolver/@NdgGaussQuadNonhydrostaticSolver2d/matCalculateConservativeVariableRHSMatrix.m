function [ termX, termY ] = matCalculateConservativeVariableRHSMatrix( obj, physClass, mesh, fphys, ftype, index)
%> @brief Function to calculate the characteristic matrix
%> @details Function to calculate the characteristic matrix in the x direction
%> @param[in] BoundaryEdge the boundary edge object
%> @param[in] InnerEdge the inner edge object
%> @param[in] Variable variable used to calculate the characteristic matrix
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[out] termX the calculated characteristic matrix in x direction
%> @param[out] termY the calculated characteristic matrix in y direction
%< Inner value and outer value of the Inner edges
[ termX, termY ] = obj.matFluxVolumeIntegral( fphys{1}(:,:,index), fphys{1}(:,:,index));

[ termX, termY ] = matEvaluateConservativeVariableTotalIntegral( obj, physClass, mesh, ftype, fphys, termX, termY, index );

termX = obj.Vq{1} * termX; termY = obj.Vq{1} * termY;

end