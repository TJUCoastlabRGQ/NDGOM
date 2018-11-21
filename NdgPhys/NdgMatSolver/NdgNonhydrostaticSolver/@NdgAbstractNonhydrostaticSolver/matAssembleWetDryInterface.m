function matAssembleWetDryInterface(obj, mesh)
%> @brief Function to assemble the wet and dry interface
%> @details
%> Function to assemble the wet and dry interface
%> @param[in] mesh The mesh object
obj.WetCellIndex = find(mesh.status == enumSWERegion.Wet);
[obj.ZeroFluxBoundary, obj.AdjacentDryCellAndFace, obj.WetDryPoint, ...
    obj.ZeroFluxBoundaryIndex, obj.WetNum] =mxAssembleWetDryInterface...
    (double(mesh.status), mesh.EToE, mesh.cell.Nface, mesh.cell.Fmask, mesh.cell.Nfp(1), mesh.cell.Np);

end