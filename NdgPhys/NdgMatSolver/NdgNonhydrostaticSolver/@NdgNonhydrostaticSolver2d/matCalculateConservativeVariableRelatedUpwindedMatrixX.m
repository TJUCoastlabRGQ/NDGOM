function [heightx, hux] = matCalculateConservativeVariableRelatedUpwindedMatrixX(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype)
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
[hm, hp] = obj.matGetFaceValue(fm(:,:,1), fp(:,:,1), ftype);
[hum, hup] = obj.matGetFaceValue(fm(:,:,2), fp(:,:,2), ftype);
[hvm, hvp] = obj.matGetFaceValue(fm(:,:,3), fp(:,:,3), ftype);
%< Inner edge contribution
fluxHmX = InnerEdge.nx.*hm;fluxHumX = InnerEdge.nx.*hum;
fluxHpX = InnerEdge.nx.*hp;fluxHupX = InnerEdge.nx.*hup;

[fluxNHx, fluxNHux] = obj.matGetUpwindedNumFluxTermX(InnerEdge, hm, hp, hum, hup, hvm, hvp);

termHx = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHmX, fluxHpX, fluxNHx );
termHux = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHumX, fluxHupX, fluxNHux );



[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys ); 
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1} );
hm = fm(:,:,1); hum = fm(:,:,2); hvm = fm(:,:,3); 
hp = fp(:,:,1); hup = fp(:,:,2); hvp = fp(:,:,3);
%< Boundary edge contribution

fluxHmX = BoundaryEdge.nx .* hm;fluxHumX = BoundaryEdge.nx .* hum;

% fluxNHx = BoundaryEdge.nx .* (hm + hp)./2;
% fluxNHux = BoundaryEdge.nx .* (hum + hup)./2;

[fluxNHx, fluxNHux] = obj.matGetUpwindedNumFluxTermX(BoundaryEdge, hm, hp, hum, hup, hvm, hvp);

termHx = - termHx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHmX, fluxNHx);
termHux = - termHux - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHumX, fluxNHux);

heightx = termHx + mesh.rx .* (mesh.cell.Dr * fphys{1}(:,:,1))...
    + mesh.sx .* (mesh.cell.Ds * fphys{1}(:,:,1));
hux = termHux + mesh.rx .* (mesh.cell.Dr * fphys{1}(:,:,2))...
    + mesh.sx .* (mesh.cell.Ds * fphys{1}(:,:,2));

end