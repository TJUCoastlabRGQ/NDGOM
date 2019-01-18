function termX = matCalculateCharacteristicMatrixX(obj, mesh,  BoundaryEdge, InnerEdge, Variable, ftype)
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
fluxMX = InnerEdge.nx.*fm;
fluxPX = InnerEdge.nx.*fp; 
termX = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMX, fluxPX);

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( Variable );        
fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, ftype, obj.EidBoundaryType);
%% test first




%< Boundary edge contribution
fluxMX = BoundaryEdge.nx.*fm;
fluxSX = BoundaryEdge.nx.*(fp + fm)./2; 
termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);

termX = termX + obj.matGetVolumeIntegralX(mesh, cell2mat(Variable));

end