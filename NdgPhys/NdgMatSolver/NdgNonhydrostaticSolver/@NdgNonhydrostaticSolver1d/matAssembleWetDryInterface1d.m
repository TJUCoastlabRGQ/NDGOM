function matAssembleWetDryInterface1d(obj, mesh)
%> @brief Function to assemble the wet and dry interface
%> @details
%> Function to assemble the wet and dry interface
%> @param[in] mesh The mesh object

[ obj.WetDryCell ] = mxAssembleWetDryInterface1d(mesh.status, mesh.EToE, ...
    mesh.InnerEdge.FToE, mesh.BoundaryEdge.FToE, mesh.BoundaryEdge.ftype, mesh.BoundaryEdge.FToF);


end