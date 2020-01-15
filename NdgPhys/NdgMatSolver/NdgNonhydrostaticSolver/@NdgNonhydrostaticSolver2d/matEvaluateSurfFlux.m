function [ fluxMPx, fluxMPy ] = matEvaluateSurfFlux( obj, status, EdgeFToE, vfmx, vfmy, Edgenx, Edgeny)
%> @brief Function to get the local flux when calculating the physical variable related first order derivative
%> @details
%> Function to get the local flux when calculating the physical variable related first order derivative,
%> we point out that this flux is set to be zero at the wet-dry interface, such that the complex partial-wet
%> cell problem can be avoided
%> @param[in] status The mesh cell status
%> @param[in] EdgeFToE The face to edge topological relations of the studied edge object
%> @param[in] vfx local or adjacent face value of the studied variable in x direction
%> @param[in] vfy local or adjacent face value of the studied variable in y direction
%> @param[in] Edgenx Projection of the unit direction vector in x direction of the studied edge object
%> @param[in] Edgeny Projection of the unit direction vector in y direction of the studied edge object
%> @param[out] fluxMPx the local or adjacent flux in x direction
%> @param[out] fluxMPy the local or adjacent flux in y direction
[ fluxMPx, fluxMPy ] = mxEvaluateSurfFlux( status, EdgeFToE, vfmx, vfmy, Edgenx, Edgeny);
end