function fluxS = matGetPrimitiveVariableInnerEdgeFlux( obj, WetDryFaceOrder, fluxS, Nfp)
%> @brief Function to calculate the numerical flux of the inner edge corresponding to the primitive variable
%> with LDG method adopted
%> @details Function to calculate the numerical flux of the inner edge corresponding to the primitive variable
%> with LDG method adopted, this version is programmed following the numerical flux for LDG method provided
%> by Peraire and Perssony (2008)
%> @param[in] WetDryFaceOrder order of wet-dry interface
%> @param[in] fluxS the numerical flux without wet-dry interface taking into consideration
%> @param[in] Nfp number of interpolation points per face
%> @param[out] fluxS the numerical flux with wet-dry interface taking into consideration
fluxS = mxGetPrimitiveVariableInnerEdgeFlux( WetDryFaceOrder, fluxS, Nfp);
end