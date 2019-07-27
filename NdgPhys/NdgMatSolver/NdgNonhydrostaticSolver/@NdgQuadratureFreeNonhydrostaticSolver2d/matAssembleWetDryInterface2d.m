function matAssembleWetDryInterface2d(obj, mesh)
%> @brief Function to assemble the wet and dry interface
%> @details
%> Function to assemble the wet and dry interface
%> @param[in] mesh The mesh object
% obj.WetCellIndex = find(mesh.status == enumSWERegion.Wet);
% [obj.ZeroFluxBoundary, obj.AdjacentDryCellAndFace, obj.WetDryPoint, ...
%     obj.ZeroFluxBoundaryIndex, obj.WetNum] =mxAssembleWetDryInterface...
%     (double(mesh.status), mesh.EToE, mesh.cell.Nface, mesh.cell.Fmask, mesh.cell.Nfp(1), mesh.cell.Np);

[ obj.WetDryCell ] = mxAssembleWetDryInterface2d(mesh.status, mesh.EToE, ...
    mesh.InnerEdge.FToE, mesh.BoundaryEdge.FToE, mesh.BoundaryEdge.ftype, mesh.BoundaryEdge.FToF);


end