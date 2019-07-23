function fluxS = matGetAuxialaryVariableInnerEdgeFlux( obj, fluxS, Um, Up, Sigmam, Sigmap, vector)
%> @brief Function to calculate the numerical flux of the inner edge corresponding to the auxialary variable
%> with LDG method adopted
%> @details Function to calculate the numerical flux of the inner edge corresponding to the auxialary variable
%> with LDG method adopted, this version is programmed following the numerical flux for LDG method provided
%> by Peraire and Perssony (2008)
%> @param[in] fluxS the numerical flux without wet-dry interface taking into consideration
%> @param[in] Um value of the local non-hydrostatic pressure
%> @param[in] Up value of the adjacent non-hydrostatic pressure
%> @param[in] Sigmam value of the local auxialary pressure
%> @param[in] Sigmap value of the adjacent auxialary pressure
%> @param[in] vector the directional vector in x or y direction
%> @param[out] fluxS the numerical flux with wet-dry interface taking into consideration

[ fluxS ] = mxGetAuxialaryVariableInnerEdgeFlux( obj.WetDryFaceOrder, fluxS, obj.Tau, ...
    obj.NonhydroFmPoint, obj.NonhydroFpPoint, Um, Up, Sigmam, Sigmap, vector);

end