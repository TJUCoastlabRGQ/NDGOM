function fluxS = matGetAuxialaryVariableBoundaryEdgeFlux( obj,  Um, Sigmam, vector)
%> @brief Function to calculate the numerical flux of the boundary edge corresponding to the auxialary variable
%> with LDG method adopted
%> @details Function to calculate the numerical flux of the boundary edge corresponding to the auxialary variable
%> with LDG method adopted, this version is programmed following the numerical flux for LDG method provided
%> by Peraire and Perssony (2008)
%> @param[in] Um value of the local non-hydrostatic pressure
%> @param[in] Sigmam value of the local auxialary pressure
%> @param[in] vector the directional vector in x or y direction
%> @param[out] fluxS the numerical flux with boundary contion taking into consideration
fluxS = mxGetAuxialaryVariableBoundaryEdgeFlux( obj.EidBoundaryType, obj.BETau, Um, Sigmam, vector);
end