function matAssembleElementBoundaryCondition(obj, mesh)
%> @brief Function to calculate the penalty parameter for the IPDG method
%> @details
%> Function to calculate the penalty parameter for the IPDG method
%> @param[in] mesh The mesh object
BoundaryEdge = mesh.BoundaryEdge;
obj.edgeType = mxAssembleElementBoundaryCondition( mesh.EToE,...
    BoundaryEdge.FToF, BoundaryEdge.FToE, int8(BoundaryEdge.ftype));
end