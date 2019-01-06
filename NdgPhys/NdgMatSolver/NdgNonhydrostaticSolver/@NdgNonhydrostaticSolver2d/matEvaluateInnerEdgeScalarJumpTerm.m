function [penaltyX, penaltyY] = matEvaluateInnerEdgeScalarJumpTerm( obj, InnerEdge, variable, eidtype)
%> @brief Function to calculate the jump value for the scalar quantity
%> @details
%> Function to calculate variable with the upwinded flux to weakly impose the continuity condition
%> @param[in] InnerEdge The inner edge object
%> @param[in] variable The variable to be calculated the jump term
%> @param[in] eidtype The boundary condition to be imposed at the wet-dry interface and all the clamped-type boundary
%> @param[out] penaltyX The penalty term in the x direction
%> @param[out] penaltyY The penalty term in the y direction
[fm, fp] = InnerEdge.matEvaluateSurfValue( variable );       
[fm, fp] = obj.matGetFaceValue(fm, fp, eidtype);

penaltyX = InnerEdge.nx .* fm - InnerEdge.nx .* fp;
penaltyY = InnerEdge.ny .* fm - InnerEdge.ny .* fp;
end