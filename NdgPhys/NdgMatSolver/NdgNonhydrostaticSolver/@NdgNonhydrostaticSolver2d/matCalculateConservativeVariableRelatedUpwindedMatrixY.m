function [heighty, huy, hvy] = matCalculateConservativeVariableRelatedUpwindedMatrixY(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype)
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
fluxHmY = InnerEdge.ny.*hm;fluxHumY = InnerEdge.ny.*hum;fluxHvmY = InnerEdge.ny.*hvm;
fluxHpY = InnerEdge.ny.*hp;fluxHupY = InnerEdge.ny.*hup;fluxHvpY = InnerEdge.ny.*hvp; 

% fluxNHy = zeros(size(hm)); fluxNHuy = zeros(size(hum)); fluxNHvy = zeros(size(hvm));
fluxNHy = InnerEdge.ny .* (hm + hp)./2; fluxNHuy = InnerEdge.ny .* (hum + hup)./2; fluxNHvy = InnerEdge.ny .* (hvm + hvp)./2;

Index = ((InnerEdge.nx .* hum + InnerEdge.ny .* hvm) > 0 & (-InnerEdge.nx .* hup - InnerEdge.ny .* hvp) < 0 );    % flow out
% fluxNHy(Index) = InnerEdge.ny(Index).* hm(Index);
fluxNHuy(Index) = InnerEdge.ny(Index).* hum(Index); fluxNHvy(Index) = InnerEdge.ny(Index).* hvm(Index);

Index = ((InnerEdge.nx .* hum + InnerEdge.ny .* hvm) < 0 & (-InnerEdge.nx .* hup - InnerEdge.ny .* hvp) > 0 );   % flow in
% fluxNHy(Index) = InnerEdge.ny(Index).* hp(Index); 
fluxNHuy(Index) = InnerEdge.ny(Index).* hup(Index); fluxNHvy(Index) = InnerEdge.ny(Index).* hvp(Index);

termHy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHmY, fluxHpY, fluxNHy );
termHuy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHumY, fluxHupY, fluxNHuy );
termHvy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHvmY, fluxHvpY, fluxNHvy );



[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys ); 
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1} );
hm = fm(:,:,1); hum = fm(:,:,2); hvm = fm(:,:,3); 
hp = fp(:,:,1); hup = fp(:,:,2); hvp = fp(:,:,3); 
%< Boundary edge contribution

fluxHmY = BoundaryEdge.ny .* hm;fluxHumY = BoundaryEdge.ny .* hum;fluxHvmY = BoundaryEdge.ny .* hvm;

% fluxNHy = zeros(size(hm)); fluxNHuy = zeros(size(hum)); fluxNHvy = zeros(size(hvm));
fluxNHy = BoundaryEdge.ny .* (hm + hp)./2; fluxNHuy = BoundaryEdge.ny .* (hum + hup)./2; fluxNHvy = BoundaryEdge.ny .* (hvm + hvp)./2;


Index = ((BoundaryEdge.nx .* hum + BoundaryEdge.ny .* hvm) > 0 & (-BoundaryEdge.nx .* hup - BoundaryEdge.ny .* hvp) < 0 );    % flow out
% fluxNHy(Index) = BoundaryEdge.ny(Index).* hm(Index); 
fluxNHuy(Index) = BoundaryEdge.ny(Index).* hum(Index); fluxNHvy(Index) = BoundaryEdge.ny(Index).* hvm(Index);

Index = ((BoundaryEdge.nx .* hum + BoundaryEdge.ny .* hvm) < 0 & (-BoundaryEdge.nx .* hup - BoundaryEdge.ny .* hvp) > 0 );   % flow in
% fluxNHy(Index) = BoundaryEdge.ny(Index).* hp(Index); 
fluxNHuy(Index) = BoundaryEdge.ny(Index).* hup(Index); fluxNHvy(Index) = BoundaryEdge.ny(Index).* hvp(Index);

termHy = - termHy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHmY, fluxNHy);
termHuy = - termHuy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHumY, fluxNHuy);
termHvy = - termHvy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHvmY, fluxNHvy);

heighty = termHy + mesh.ry .* (mesh.cell.Dr * fphys{1}(:,:,1))...
    + mesh.sy .* (mesh.cell.Ds * fphys{1}(:,:,1));
huy = termHuy + mesh.ry .* (mesh.cell.Dr * fphys{1}(:,:,2))...
    + mesh.sy .* (mesh.cell.Ds * fphys{1}(:,:,2));
hvy = termHvy + mesh.ry .* (mesh.cell.Dr * fphys{1}(:,:,3))...
    + mesh.sy .* (mesh.cell.Ds * fphys{1}(:,:,3));

end