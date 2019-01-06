function [termx, termy] = matCalculateLDGTerm( obj, mesh, BoundaryEdge, InnerEdge, Variable, qx, qy, eidtype)

tau = 1.0;
Beta = [1, 1];

%< Inner edge contribution
[penaltyX, penaltyY] = obj.matEvaluateInnerEdgeScalarJumpTerm( InnerEdge, Variable, eidtype);  %Penalty term for u
[ penaltyTerm ] =  obj.matEvaluateInnerEdgeVectorJumpTerm( InnerEdge, qx, qy, eidtype);


[fmx, fpx] = InnerEdge.matEvaluateSurfValue( qx );  [fmy, fpy] = InnerEdge.matEvaluateSurfValue( qy );        
[fmx, fpx] = obj.matGetFaceValue(fmx, fpx, eidtype); [fmy, fpy] = obj.matGetFaceValue(fmy, fpy, eidtype);


fluxMX = InnerEdge.nx.*fmx; fluxMY = InnerEdge.ny.*fmy;
fluxPX = InnerEdge.nx.*fpx; fluxPY = InnerEdge.ny.*fpy; 

fluxSx = InnerEdge.nx .* ( (fmx + fpx)./2 + Beta(1) .* penaltyTerm  - tau .* penaltyX);
fluxSy = InnerEdge.ny .* ( (fmy + fpy)./2 + Beta(2) .* penaltyTerm  - tau .* penaltyY);

termx = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMX, fluxPX, fluxSx );
termy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMY, fluxPY, fluxSy );

%< boundary edge contribution
[penaltyX, penaltyY] = obj.matEvaluateBoundaryEdgeScalarJumpTerm( BoundaryEdge, Variable, eidtype);  %Penalty term for u
[ penaltyTerm ] =  obj.matEvaluateBoundaryEdgeVectorJumpTerm( BoundaryEdge, qx, qy, eidtype);


[fmx, fpx] = BoundaryEdge.matEvaluateSurfValue( qx ); [fmy, fpy] = BoundaryEdge.matEvaluateSurfValue( qy );         
fpx = obj.matImposeNonhydroRelatedBoundaryCondition(fmx, fpx, eidtype, obj.EidBoundaryType);
fpy = obj.matImposeNonhydroRelatedBoundaryCondition(fmy, fpy, eidtype, obj.EidBoundaryType);
%< Boundary edge contribution
fluxMX = BoundaryEdge.nx.*fmx; fluxMY = BoundaryEdge.ny.*fmy;

fluxSX = BoundaryEdge.nx .* ( (fmx + fpx)./2 + Beta(1) .* penaltyTerm  - tau .* penaltyX);
fluxSY = BoundaryEdge.ny .* ( (fmy + fpy)./2 + Beta(2) .* penaltyTerm  - tau .* penaltyY);

termx = - termx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
termy = - termy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);

termx = termx + mesh.rx .* (mesh.cell.Dr * cell2mat(qx))...
    + mesh.sx .* (mesh.cell.Ds * cell2mat(qx));
termy = termy + mesh.ry .* (mesh.cell.Dr * cell2mat(qy))...
    + mesh.sy .* (mesh.cell.Ds * cell2mat(qy));


end