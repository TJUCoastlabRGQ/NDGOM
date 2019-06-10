function [Um, Up, fluxS] = matGetWetDryFaceVariableAndFlux( obj, Um, Up, fluxS)
%> @brief Function to get the variable value and the numerical flux at the wet-dry interface
%> @details Function to get the variable value and the numerical flux at the wet-dry interface
%> @param[in] Um value of the local cell
%> @param[in] Up value of the adjacent cell
%> @param[in] fluxS the flux without wet-dry interface considered
%> @param[out] Um value of the local cell with wet-dry interface considered
%> @param[out] Up value of the adjacent cell with wet-dry interface considered
%> @param[in] fluxS the flux with wet-dry interface considered
[Um, Up, fluxS] = mxGetWetDryFaceVariableAndFlux( obj.WetDryFaceOrder, obj.NonhydroFmPoint, obj.NonhydroFpPoint, Um, Up, fluxS );
end