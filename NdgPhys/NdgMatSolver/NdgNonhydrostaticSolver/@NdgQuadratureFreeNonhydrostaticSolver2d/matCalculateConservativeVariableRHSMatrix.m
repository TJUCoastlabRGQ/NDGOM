function [ termX, termY ] = matCalculateConservativeVariableRHSMatrix( obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype, index)
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
fluxMX = InnerEdge.nx.*fm;  fluxMY = InnerEdge.ny.*fm;
fluxPX = InnerEdge.nx.*fp;  fluxPY = InnerEdge.ny.*fp;
termX = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMX, fluxPX);
termY = InnerEdge.matEvaluateStrongFormEdgeCentralRHS(fluxMY, fluxPY);


[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys ); 
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1} );

% for i = 1:numel(BoundaryEdge.ftype)
%     if BoundaryEdge.ftype(i) == enumBoundaryCondition.ClampedVel
%         fp(:,i,:) = fm(:,i,:);
%     end
% end

%< Boundary edge contribution
fluxMX = BoundaryEdge.nx.*fm(:,:,index); fluxMY = BoundaryEdge.ny.*fm(:,:,index);
fluxSX = BoundaryEdge.nx.*(fp(:,:,index) + fm(:,:,index))./2; fluxSY = BoundaryEdge.ny.*(fp(:,:,index) + fm(:,:,index))./2;
termX = - termX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
termY = - termY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);

[VolumeIntegralX, VolumeIntegralY] = obj.matVolumeIntegral( mesh, fphys{1}(:,:,index), fphys{1}(:,:,index));

termX = termX + VolumeIntegralX;
termY = termY + VolumeIntegralY;
end