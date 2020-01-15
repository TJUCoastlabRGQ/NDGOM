function faceflag = matGetWetDryFace( obj, status, InnerEdgeFToE)
%> @brief Function to get the inner edge that cells located besides are wet and non-wet cell
%> @details Function to get the inner edge that cells located besides are wet and non-wet cell
%> @param[in] status Status of the mesh cell
%> @param[in] InnerEdgeFToE The face to element topological relation of the inner edge
%> @param[out] faceflag the face status flag, 1 for wet-dry interface, 0 for wet-wet interface
faceflag = mxGetWetDryFace( status, InnerEdgeFToE);
end