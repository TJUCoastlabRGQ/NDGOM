function fluxS = matGetPrimitiveVariableBoundaryEdgeFlux( obj, vector, Um )
%> @brief Function to calculate the numerical flux of the boundary edge corresponding to the primitive variable
%> with LDG method adopted
%> @details Function to calculate the numerical flux of the boundary edge corresponding to the primitive variable
%> with LDG method adopted, this version is programmed following the numerical flux for LDG method provided
%> by Peraire and Perssony (2008)
%> @param[in] vector the directional vector in x or y direction
%> @param[in] Um value of the local non-hydrostatic pressure
%> @param[out] fluxS the numerical flux in direction x or y
fluxS = mxGetPrimitiveVariableBoundaryEdgeFlux( obj.EidBoundaryType, vector, Um );
end