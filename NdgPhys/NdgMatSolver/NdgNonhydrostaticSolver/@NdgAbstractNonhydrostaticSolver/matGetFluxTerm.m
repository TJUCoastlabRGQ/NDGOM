function FluxTerm = matGetFluxTerm(obj, mesh, FluxTerm)
%< @brief This function is used to get the flux term at the wet dry interface
%> @details
%> Function to get the flux term at the wet dry interface
%> @param[in] physClass The hydrostatic solver
%> @param[in] FluxTerm The flux term needs to be corrected
%> @param[in] mesh The mesh object
%> @param[out] FluxTerm The flux term corrected
FluxTerm = mxGetFluxTerm(obj.ZeroFluxBoundary, obj.ZeroFluxBoundaryIndex, ...
    obj.AdjacentDryCellAndFace, mesh.cell.Nfp(1), FluxTerm, mesh.cell.TNfp);

end