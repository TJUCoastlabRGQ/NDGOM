function [ fluxMPx ] = matEvaluateSurfFlux( obj, status, EdgeFToE, vfmx, Edgenx )
%> @brief Function to get the local flux when calculating the physical variable related first order derivative
%> @details
%> Function to get the local flux when calculating the physical variable related first order derivative,
%> we point out that this flux is set to be zero at the wet-dry interface, such that the complex partial-wet
%> cell problem can be avoided
%> @param[in] status The mesh cell status
%> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
%> @param[in] vfmx local or adjacent face value of the studied variable in x direction
%> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
%> @param[out] fluxMPx the local or adjacent flux in x direction
[ fluxMPx ] = mxEvaluateSurfFlux( status, EdgeFToE, vfmx, Edgenx );
end