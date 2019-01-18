function RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, physClass)

InnerEdge = physClass.meshUnion(1).InnerEdge;
BoundaryEdge = physClass.meshUnion(1).BoundaryEdge;
fhx = obj.matCalculateConservativeVariableRelatedMatrixX(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 1);
fhy = obj.matCalculateConservativeVariableRelatedMatrixY(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 1);
fhux = obj.matCalculateConservativeVariableRelatedMatrixX(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 2);
fhvy = obj.matCalculateConservativeVariableRelatedMatrixY(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 3);

% [fhx, fhux] = obj.matCalculateConservativeVariableRelatedUpwindedMatrixX( physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero);
% [fhy, fhvy] = obj.matCalculateConservativeVariableRelatedUpwindedMatrixY( physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero);

% fhux = obj.matCalculateUpwindedConservativeVariableRelatedMatrixX(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 2);
% fhvy = obj.matCalculateUpwindedConservativeVariableRelatedMatrixY(physClass, BoundaryEdge, InnerEdge, fphys, enumNonhydroBoundaryCondition.Zero, 3);
% RHS = -2 * ((fphys{1}(:,:,6) - fphys{1}(:,:,2) .* obj.bx  - fphys{1}(:,:,3) .* obj.by)) + ...
%      ((fhx.* fphys{1}(:,:,2))) + ((fhy.* fphys{1}(:,:,3)))  - ((fhux.* fphys{1}(:,:,1))) - ...
%      ((fhvy.* fphys{1}(:,:,1)));

RHS = -2 * fphys{1}(:,:,6) + 2 * fphys{1}(:,:,2) .* obj.bx  + 2 * fphys{1}(:,:,3) .* obj.by + ...
  fhx.* fphys{1}(:,:,2) + fhy.* fphys{1}(:,:,3)  - fhux.* fphys{1}(:,:,1) - fhvy.* fphys{1}(:,:,1);
 
% [BFIx, BFIy] = getBoundaryFaceIntegral(physClass, BoundaryEdge);
% [BGFIx, BGFIy] = resembleAndProductBoundaryFaceIntegral(physClass, BFIx, BFIy, BoundaryEdge);
 
 %% New version, 
%  [~, Hx, Hy] = getDepthRelatedTerm(fphys, physClass);
 
%  RHS = RHS + obj.dt*H.*(mesh.rx.*(mesh.cell.Dr*BGFIx)+mesh.sx.*(mesh.cell.Ds*BGFIx))+...
%      obj.dt*H.*(mesh.ry.*(mesh.cell.Dr*BGFIy)+mesh.sy.*(mesh.cell.Ds*BGFIy))-...
%      obj.dt*Phx.*BGFIx - obj.dt*Phy.*BGFIy;
% 
%  RHS = RHS + obj.dt*H.*(mesh.rx.*(mesh.cell.Dr*BGFIx)+mesh.sx.*(mesh.cell.Ds*BGFIx))+...
%      obj.dt*H.*(mesh.ry.*(mesh.cell.Dr*BGFIy)+mesh.sy.*(mesh.cell.Ds*BGFIy))-...
%      obj.dt*Hx.*BGFIx - obj.dt*Hy.*BGFIy;

%  RHS = RHS - obj.dt*Hx.*BGFIx - obj.dt*Hy.*BGFIy;
end


function [H, Hx, Hy] = getDepthRelatedTerm(fphys, physClass)

H = fphys{1}(:,:,1);

mesh = physClass.meshUnion(1);
InnerEdge = physClass.meshUnion(1).InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue(fphys);
fluxHxm = InnerEdge.nx .* fm(:,:,1);  fluxHym = InnerEdge.ny .* fm(:,:,1);
fluxHxp = InnerEdge.nx .* fp(:,:,1); fluxHyp =  InnerEdge.ny .* fp(:,:,1);
HrhsX = InnerEdge.matEvaluateStrongFormEdgeCentralRHS( fluxHxm, fluxHxp);
HrhsY = InnerEdge.matEvaluateStrongFormEdgeCentralRHS( fluxHym, fluxHyp);

BoundaryEdge = physClass.meshUnion(1).BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue(fphys);
[fm, fp] = physClass.matImposeBoundaryCondition(BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1});
fluxMX = BoundaryEdge.nx.*fm(:,:,1);  fluxMY = BoundaryEdge.ny.*fm(:,:,1);
fluxSX = BoundaryEdge.nx.*(fp(:,:,1) + fm(:,:,1))./2;   fluxSY = BoundaryEdge.ny.*(fp(:,:,1) + fm(:,:,1))./2; 
HrhsX = - HrhsX - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
HrhsY = - HrhsY - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);

Hx = HrhsX + mesh.rx .* (mesh.cell.Dr * H)...
    + mesh.sx .* (mesh.cell.Ds * H);
Hy = HrhsY + mesh.ry .* (mesh.cell.Dr * H)...
    + mesh.sy .* (mesh.cell.Ds * H);
end



function [XBoundaryFaceIntegra, YBoundaryFaceIntegra] = getBoundaryFaceIntegral(physClass, BoundaryEdge)
Ext = physClass.fext{1}(:,:,6);
XBoundaryFaceIntegra = zeros(BoundaryEdge.Nfp, BoundaryEdge.Ne);
YBoundaryFaceIntegra = zeros(BoundaryEdge.Nfp, BoundaryEdge.Ne);
XDirecExt = Ext .* BoundaryEdge.nx;
YDirecExt = Ext .* BoundaryEdge.ny; % with boundary normal direction considered
for i = 1:BoundaryEdge.Ne
    XBoundaryFaceIntegra(:,i) = BoundaryEdge.M*diag(BoundaryEdge.Js(:,i))*XDirecExt(:,i);
    YBoundaryFaceIntegra(:,i) = BoundaryEdge.M*diag(BoundaryEdge.Js(:,i))*YDirecExt(:,i);
end
end


function [BGFIx, BGFIy] = resembleAndProductBoundaryFaceIntegral(physClass, BFIx, BFIy, BoundaryEdge)
invM = physClass.meshUnion(1).cell.invM;
BGFIx = zeros(size(physClass.fphys{1}(:,:,1)));
BGFIy = zeros(size(physClass.fphys{1}(:,:,1)));
for i = 1:BoundaryEdge.Ne
    ele = BoundaryEdge.FToE(1,i);
    for j = 1:BoundaryEdge.Nfp
        BGFIx(BoundaryEdge.FToN1(j,ele),ele) =  BGFIx(BoundaryEdge.FToN1(j,ele),ele) + BFIx(j,i);
        BGFIy(BoundaryEdge.FToN1(j,ele),ele) =  BGFIy(BoundaryEdge.FToN1(j,ele),ele) + BFIy(j,i);
    end
end
BGFIx = 1./physClass.meshUnion(1).J .* (invM * BGFIx);
BGFIy = 1./physClass.meshUnion(1).J .* (invM * BGFIy);
end