function [ penaltyTerm ] = matEvaluateBoundaryEdgeVectorJumpTerm( obj, BoundaryEdge, variablex, variabley, eidtype)
%> @brief Function to calculate the jump value for the scalar quantity
%> @details
%> Function to calculate variable with the upwinded flux to weakly impose the continuity condition
%> @param[in] InnerEdge The inner edge object
%> @param[in] variable The variable to be calculated the jump term
%> @param[in] eidtype The boundary condition to be imposed at the wet-dry interface and all the clamped-type boundary
%> @param[out] penaltyX The penalty term in the x direction
%> @param[out] penaltyY The penalty term in the y direction
[fmx, fpx] = BoundaryEdge.matEvaluateSurfValue( variablex );  [fmy, fpy] = BoundaryEdge.matEvaluateSurfValue( variabley );        
fpx = obj.matImposeNonhydroRelatedBoundaryCondition(fmx, fpx, eidtype, obj.EidBoundaryType);
fpy = obj.matImposeNonhydroRelatedBoundaryCondition(fmy, fpy, eidtype, obj.EidBoundaryType);

penaltyTerm = BoundaryEdge.nx .* fmx - BoundaryEdge.nx .* fpx + BoundaryEdge.ny .* fmy - BoundaryEdge.ny .* fpy;
end