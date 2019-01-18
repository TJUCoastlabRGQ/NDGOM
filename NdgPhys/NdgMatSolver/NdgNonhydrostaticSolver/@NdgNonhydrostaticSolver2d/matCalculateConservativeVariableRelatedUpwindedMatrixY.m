function [heighty, hvy] = matCalculateConservativeVariableRelatedUpwindedMatrixY(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype)
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
fluxHmY = InnerEdge.ny.*hm;fluxHvmY = InnerEdge.ny.*hvm;
fluxHpY = InnerEdge.ny.*hp;fluxHvpY = InnerEdge.ny.*hvp; 

[fluxNHy, fluxNHvy] = obj.matGetUpwindedNumFluxTermY(InnerEdge, hm, hp, hum, hup, hvm, hvp);

termHy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHmY, fluxHpY, fluxNHy );
termHvy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHvmY, fluxHvpY, fluxNHvy );

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys ); 
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1} );
hm = fm(:,:,1); hum = fm(:,:,2);  hvm = fm(:,:,3); 
hp = fp(:,:,1); hup = fp(:,:,2);  hvp = fp(:,:,3); 
%< Boundary edge contribution

fluxHmY = BoundaryEdge.ny .* hm;fluxHvmY = BoundaryEdge.ny .* hvm;
% fluxNHy = BoundaryEdge.ny .* (hm + hp)./2;
% fluxNHvy = BoundaryEdge.ny .* (hvm + hvp)./2;

[fluxNHy, fluxNHvy] = obj.matGetUpwindedNumFluxTermY(BoundaryEdge, hm, hp, hum, hup, hvm, hvp);


termHy = - termHy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHmY, fluxNHy);
termHvy = - termHvy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHvmY, fluxNHvy);

heighty = termHy + mesh.ry .* (mesh.cell.Dr * fphys{1}(:,:,1))...
    + mesh.sy .* (mesh.cell.Ds * fphys{1}(:,:,1));
hvy = termHvy + mesh.ry .* (mesh.cell.Dr * fphys{1}(:,:,3))...
    + mesh.sy .* (mesh.cell.Ds * fphys{1}(:,:,3));

end