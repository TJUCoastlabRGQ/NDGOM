function [heightx, hux, hvx] = matCalculateConservativeVariableRelatedUpwindedMatrixX(obj, physClass, BoundaryEdge, InnerEdge, fphys, ftype)
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
fluxHmX = InnerEdge.nx.*hm;fluxHumX = InnerEdge.nx.*hum;fluxHvmX = InnerEdge.nx.*hvm;
fluxHpX = InnerEdge.nx.*hp;fluxHupX = InnerEdge.nx.*hup;fluxHvpX = InnerEdge.nx.*hvp; 

hu_norm = InnerEdge.ny .* hum;
sign_hu_norm = sign(hu_norm);
fluxNHx = InnerEdge.nx .* (hm.*(sign_hu_norm+1)*0.5 + hp.*(1 - sign_hu_norm)*0.5) .* sign_hu_norm;
fluxNHux = InnerEdge.nx .* (hum.*(sign_hu_norm+1)*0.5 + hup.*(1 - sign_hu_norm)*0.5) .* sign_hu_norm;
fluxNHvx = InnerEdge.nx .* (hvm.*(sign_hu_norm+1)*0.5 + hvp.*(1 - sign_hu_norm)*0.5) .* sign_hu_norm;

% % fluxNHx = zeros(size(hm)); fluxNHux = zeros(size(hum)); fluxNHvx = zeros(size(hvm));
% fluxNHx = InnerEdge.nx .* (hm + hp)./2; fluxNHux = InnerEdge.nx .* (hum + hup)./2; fluxNHvx = InnerEdge.nx .* (hvm + hvp)./2;
% 
% Index = ((InnerEdge.nx .* hum + InnerEdge.ny .* hvm) > 0 & (-InnerEdge.nx .* hup - InnerEdge.ny .* hvp) < 0 );    % flow out
% % fluxNHx(Index) = InnerEdge.nx(Index).* hm(Index);
% fluxNHux(Index) = InnerEdge.nx(Index).* hum(Index); fluxNHvx(Index) = InnerEdge.nx(Index).* hvm(Index);
% 
% Index = ((InnerEdge.nx .* hum + InnerEdge.ny .* hvm) < 0 & (-InnerEdge.nx .* hup - InnerEdge.ny .* hvp) > 0 );   % flow in
% % fluxNHx(Index) = InnerEdge.nx(Index).* hp(Index);
% fluxNHux(Index) = InnerEdge.nx(Index).* hup(Index); fluxNHvx(Index) = InnerEdge.nx(Index).* hvp(Index);


termHx = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHmX, fluxHpX, fluxNHx );
termHux = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHumX, fluxHupX, fluxNHux );
termHvx = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxHvmX, fluxHvpX, fluxNHvx );


[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys ); 
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1} );
hm = fm(:,:,1); hum = fm(:,:,2); hvm = fm(:,:,3); 
hp = fp(:,:,1); hup = fp(:,:,2); hvp = fp(:,:,3); 
%< Boundary edge contribution

fluxHmX = BoundaryEdge.nx .* hm;fluxHumX = BoundaryEdge.nx .* hum;fluxHvmX = BoundaryEdge.nx .* hvm;

% % fluxNHx = zeros(size(hm)); fluxNHux = zeros(size(hum)); fluxNHvx = zeros(size(hvm));
% fluxNHx = BoundaryEdge.nx .* (hm + hp)./2; fluxNHux = BoundaryEdge.nx .* (hum + hup)./2; fluxNHvx = BoundaryEdge.nx .* (hvm + hvp)./2;
% 
% Index = ((BoundaryEdge.nx .* hum + BoundaryEdge.ny .* hvm) > 0 & (-BoundaryEdge.nx .* hup - BoundaryEdge.ny .* hvp) < 0 );    % flow out
% % fluxNHx(Index) = BoundaryEdge.nx(Index).* hm(Index);
% fluxNHux(Index) = BoundaryEdge.nx(Index).* hum(Index); fluxNHvx(Index) = BoundaryEdge.nx(Index).* hvm(Index);
% 
% Index = ((BoundaryEdge.nx .* hum + BoundaryEdge.ny .* hvm) < 0 & (-BoundaryEdge.nx .* hup - BoundaryEdge.ny .* hvp) > 0 );   % flow in
% % fluxNHx(Index) = BoundaryEdge.nx(Index).* hp(Index); 
% fluxNHux(Index) = BoundaryEdge.nx(Index).* hup(Index); fluxNHvx(Index) = BoundaryEdge.nx(Index).* hvp(Index);

hu_norm = BoundaryEdge.ny .* hum;
sign_hu_norm = sign(hu_norm);
fluxNHx = BoundaryEdge.nx .* (hm.*(sign_hu_norm+1)*0.5 + hp.*(1 - sign_hu_norm)*0.5) .* sign_hu_norm;
fluxNHux = BoundaryEdge.nx .* (hum.*(sign_hu_norm+1)*0.5 + hup.*(1 - sign_hu_norm)*0.5) .* sign_hu_norm;
fluxNHvx = BoundaryEdge.nx .* (hvm.*(sign_hu_norm+1)*0.5 + hvp.*(1 - sign_hu_norm)*0.5) .* sign_hu_norm;

termHx = - termHx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHmX, fluxNHx);
termHux = - termHux - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHumX, fluxNHux);
termHvx = - termHvx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxHvmX, fluxNHvx);

heightx = termHx + mesh.rx .* (mesh.cell.Dr * fphys{1}(:,:,1))...
    + mesh.sx .* (mesh.cell.Ds * fphys{1}(:,:,1));
hux = termHux + mesh.rx .* (mesh.cell.Dr * fphys{1}(:,:,2))...
    + mesh.sx .* (mesh.cell.Ds * fphys{1}(:,:,2));
hvx = termHvx + mesh.rx .* (mesh.cell.Dr * fphys{1}(:,:,3))...
    + mesh.sx .* (mesh.cell.Ds * fphys{1}(:,:,3));

end