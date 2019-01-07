function tempq2x = matCalculatePenaltyCharacteristicMatrixX(obj, mesh, BoundaryEdge, InnerEdge, gmat, tempqx, EidBoundaryType)
tau = 1.0;
%% Inner Edge part
[um, up] = InnerEdge.matEvaluateSurfValue(gmat);
[um, up] = obj.matGetFaceValue(um, up, enumNonhydroBoundaryCondition.Zero);
JumpU = InnerEdge.nx .* (um - up) + InnerEdge.ny .* (um - up);   % Calculate the jump term finished

[qm, qp] = InnerEdge.matEvaluateSurfValue(tempqx);
[qm, qp] = obj.matGetFaceValue(qm, qp, enumNonhydroBoundaryCondition.ZeroGrad);
fluxq = (qm + qp)./2.*InnerEdge.nx - tau * JumpU.*InnerEdge.nx; %Calculate the numerical flux finished

fluxM = InnerEdge.nx .* qm;   
fluxP = InnerEdge.nx .* qp;

tempq2x = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxM, fluxP, fluxq); %Inner Edge Integration finished

%% Boundary Edge part
[um, up] = BoundaryEdge.matEvaluateSurfValue(gmat);
up = obj.matImposeNonhydroRelatedBoundaryCondition(um, up, enumNonhydroBoundaryCondition.Zero, EidBoundaryType);
JumpU = BoundaryEdge.nx .* (um - up) + BoundaryEdge.ny .* (um - up);   % Calculate the jump term finished


[qm, qp] = BoundaryEdge.matEvaluateSurfValue(tempqx);
qp = obj.matImposeNonhydroRelatedBoundaryCondition(qm, qp, enumNonhydroBoundaryCondition.ZeroGrad, EidBoundaryType);
fluxq = (qm + qp)./2.*BoundaryEdge.nx - tau * JumpU.*BoundaryEdge.nx; %Calculate the numerical flux finished

fluxM = BoundaryEdge.nx .* qm; 

tempq2x = -tempq2x - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxM, fluxq);

tempq2x = tempq2x + mesh.rx .* (mesh.cell.Dr * cell2mat(tempqx))...
    + mesh.sx .* (mesh.cell.Ds * cell2mat(tempqx));
end