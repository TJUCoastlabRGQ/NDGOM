function matCalculatePenaltyParameter( obj, mesh )
%> @brief Function to calculate the penalty parameter for the IPDG method
%> @details
%> Function to calculate the penalty parameter for the IPDG method
%> @param[in] mesh The mesh object
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
obj.Tau = 2 * mxCalculatePenaltyParameter(mesh.LAV, InnerEdge.FToE, InnerEdge.LAV,...
    BoundaryEdge.FToE, BoundaryEdge.LAV, int8(mesh.type), mesh.cell.N);
end