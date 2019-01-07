function tempq2y = matCalculatePenaltyCharacteristicMatrixY(obj, mesh, BoundaryEdge, InnerEdge, gmat, tempqy, EidBoundaryType )
tau = 1.0;
%% Inner Edge part
[um, up] = InnerEdge.matEvaluateSurfValue(gmat);
[um, up] = obj.matGetFaceValue(um, up, enumNonhydroBoundaryCondition.Zero);
JumpU = InnerEdge.nx .* (um - up) + InnerEdge.ny .* (um - up);   % Calculate the jump term finished

[qm, qp] = InnerEdge.matEvaluateSurfValue(tempqy);
[qm, qp] = obj.matGetFaceValue(qm, qp, enumNonhydroBoundaryCondition.ZeroGrad);
fluxq = (qm + qp)./2.*InnerEdge.ny - tau * JumpU.*InnerEdge.ny; %Calculate the numerical flux finished

fluxM = InnerEdge.ny .* qm;   
fluxP = InnerEdge.ny .* qp;

tempq2y = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxM, fluxP, fluxq); %Inner Edge Integration finished

%% Boundary Edge part
[um, up] = BoundaryEdge.matEvaluateSurfValue(gmat);
up = obj.matImposeNonhydroRelatedBoundaryCondition(um, up, enumNonhydroBoundaryCondition.Zero, EidBoundaryType);
JumpU = BoundaryEdge.nx .* (um - up) + BoundaryEdge.ny .* (um - up);   % Calculate the jump term finished


[qm, qp] = BoundaryEdge.matEvaluateSurfValue(tempqy);
qp = obj.matImposeNonhydroRelatedBoundaryCondition(qm, qp, enumNonhydroBoundaryCondition.ZeroGrad, EidBoundaryType);
fluxq = (qm + qp)./2.*BoundaryEdge.ny - tau * JumpU.*BoundaryEdge.ny; %Calculate the numerical flux finished

fluxM = BoundaryEdge.ny .* qm; 

tempq2y = -tempq2y - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxM, fluxq);

tempq2y = tempq2y + mesh.ry .* (mesh.cell.Dr * cell2mat(tempqy))...
    + mesh.sy .* (mesh.cell.Ds * cell2mat(tempqy));
end