function termY = matCalculateCharacteristicMatrixY(obj, mesh,  BoundaryEdge, InnerEdge, Variable, ftype)
%> @brief Function to calculate the characteristic matrix
%> @details Function to calculate the characteristic matrix in the x direction
%> @param[in] BoundaryEdge the boundary edge object
%> @param[in] InnerEdge the inner edge object
%> @param[in] Variable variable used to calculate the characteristic matrix
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[out] termX the calculated characteristic matrix in x direction
%> @param[out] termY the calculated characteristic matrix in y direction
%< Inner value and outer value of the Inner edges
[fm, fp] = InnerEdge.matEvaluateSurfValue( Variable );       
[fm, fp] = obj.matGetFaceValue(fm, fp, ftype);
%< Inner edge contribution
fluxMY = InnerEdge.ny.*fm;
fluxPY = InnerEdge.ny.*fp;
termY = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMY, fluxPY);

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( Variable );        
fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, ftype, obj.EidBoundaryType);

%% test first


%< Boundary edge contribution
fluxMY = BoundaryEdge.ny.*fm; 
fluxSY = BoundaryEdge.ny.*(fp + fm)./2;
termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);


termY = termY + obj.matGetVolumeIntegralY(mesh, cell2mat(Variable));
end