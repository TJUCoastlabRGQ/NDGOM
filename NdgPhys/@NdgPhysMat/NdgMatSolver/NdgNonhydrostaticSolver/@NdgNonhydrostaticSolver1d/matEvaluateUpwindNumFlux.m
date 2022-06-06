function  [ fluxX ] = matEvaluateUpwindNumFlux( obj, status, EdgeFToE, HUm, HUp, vfmx, vfpx, Edgenx )
%> @brief Function to get the upwinded numerical flux when calculating the physical variable
%> related first order derivative in upwind manner
%> @details
%> Function to get the upwinded numerical flux when calculating the physical variable related first
%> order derivative in upwind manner, we point out that this flux is set to be zero at the wet-dry
%> interface, such that the complex partial- wet cell problem can be avoided
%> @param[in] status The mesh cell status
%> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
%> @param[in] HUm The HU field of the local face
%> @param[in] HUp The HU field of the adjacent face
%> @param[in] vfmx local face value of the studied variable in x direction
%> @param[in] vfpx adjacent face value of the studied variable in x direction
%> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
%> @param[out] fluxX the upwind numerical flux in x direction
[ fluxX ] = mxEvaluateUpwindNumFlux( status, EdgeFToE, HUm, HUp, vfmx, vfpx, Edgenx );
end