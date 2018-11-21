function termY = matCalculateConservativeVariableRelatedMatrixY(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype, index)
%> @brief Function to calculate the characteristic matrix
%> @details Function to calculate the characteristic matrix in the x direction
%> @param[in] BoundaryEdge the boundary edge object
%> @param[in] InnerEdge the inner edge object
%> @param[in] Variable variable used to calculate the characteristic matrix
%> @param[in] ftype enumeration type used to impose the non-hydro static relalated boundary condition at the wet dry interface
%> @param[out] termX the calculated characteristic matrix in x direction
%> @param[out] termY the calculated characteristic matrix in y direction
%< Inner value and outer value of the Inner edges
mesh = physClass.meshUnion(1);
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );       
[fm, fp] = obj.matGetFaceValue(fm(:,:,index), fp(:,:,index), ftype);
%< Inner edge contribution
fluxMY = InnerEdge.ny.*fm;
fluxPY = InnerEdge.ny.*fp; 
termY = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMY, fluxPY);


[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys ); 
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1} );

%< Boundary edge contribution
fluxMY = BoundaryEdge.ny.*fm(:,:,index);
fluxSX = BoundaryEdge.ny.*(fp(:,:,index) + fm(:,:,index))./2; 
termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSX);

termY = termY + mesh.rx .* (mesh.cell.Dr * fphys{1}(:,:,index))...
    + mesh.sx .* (mesh.cell.Ds * fphys{1}(:,:,index));

end