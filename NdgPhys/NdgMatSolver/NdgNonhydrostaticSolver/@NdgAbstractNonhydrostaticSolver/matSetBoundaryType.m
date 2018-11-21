function matSetBoundaryType(obj,mesh)
%> @brief Function to set the cell boundary type
%> @details
%> Function to  set the cell boundary type to impose the nonhydrostatic related boundary condition
%> @param[in] mesh The mesh object
obj.EidBoundaryType = mxSetBoundaryType(int8(mesh.BoundaryEdge.ftype),...
   mesh.BoundaryEdge.FToN1);

end