function [ penaltyTerm ] = matEvaluateInnerEdgeVectorJumpTerm( obj, InnerEdge, variablex, variabley, eidtype)
%> @brief Function to calculate the jump value for the scalar quantity
%> @details
%> Function to calculate variable with the upwinded flux to weakly impose the continuity condition
%> @param[in] InnerEdge The inner edge object
%> @param[in] variable The variable to be calculated the jump term
%> @param[in] eidtype The boundary condition to be imposed at the wet-dry interface and all the clamped-type boundary
%> @param[out] penaltyX The penalty term in the x direction
%> @param[out] penaltyY The penalty term in the y direction
[fmx, fpx] = InnerEdge.matEvaluateSurfValue( variablex ); [fmy, fpy] = InnerEdge.matEvaluateSurfValue( variabley );         
[fmx, fpx] = obj.matGetFaceValue(fmx, fpx, eidtype); [fmy, fpy] = obj.matGetFaceValue(fmy, fpy, eidtype);

penaltyTerm = InnerEdge.nx .* fmx - InnerEdge.nx .* fpx + InnerEdge.ny .* fmy - InnerEdge.ny .* fpy;
end