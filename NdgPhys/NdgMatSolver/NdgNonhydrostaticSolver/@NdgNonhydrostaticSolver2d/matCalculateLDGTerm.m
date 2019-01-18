function [termx, termy] = matCalculateLDGTerm( obj, mesh, BoundaryEdge, InnerEdge, Variable, qx, qy)

tau = 1;
Beta = [1, 1];

%< Inner edge contribution
[penaltyX, penaltyY] = obj.matEvaluateInnerEdgeScalarJumpTerm( InnerEdge, Variable, enumNonhydroBoundaryCondition.Zero);  %Penalty term for u
[ penaltyTerm ] =  obj.matEvaluateInnerEdgeVectorJumpTerm( InnerEdge, qx, qy, enumNonhydroBoundaryCondition.ZeroGrad);


[fmx, fpx] = InnerEdge.matEvaluateSurfValue( qx );  [fmy, fpy] = InnerEdge.matEvaluateSurfValue( qy );        
[fmx, fpx] = obj.matGetFaceValue(fmx, fpx, enumNonhydroBoundaryCondition.ZeroGrad); 
[fmy, fpy] = obj.matGetFaceValue(fmy, fpy, enumNonhydroBoundaryCondition.ZeroGrad);


fluxMX = InnerEdge.nx.*fmx; fluxMY = InnerEdge.ny.*fmy;
fluxPX = InnerEdge.nx.*fpx; fluxPY = InnerEdge.ny.*fpy; 

fluxSX = InnerEdge.nx .* ( (fmx + fpx)./2 + Beta(1) .* penaltyTerm  - tau .* penaltyX);
fluxSY = InnerEdge.ny .* ( (fmy + fpy)./2 + Beta(2) .* penaltyTerm  - tau .* penaltyY);

termx = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMX, fluxPX, fluxSX );
termy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMY, fluxPY, fluxSY );

%< boundary edge contribution
[penaltyX, penaltyY] = obj.matEvaluateBoundaryEdgeScalarJumpTerm( BoundaryEdge, Variable, enumNonhydroBoundaryCondition.Zero);  %Penalty term for u
[ penaltyTerm ] =  obj.matEvaluateBoundaryEdgeVectorJumpTerm( BoundaryEdge, qx, qy, enumNonhydroBoundaryCondition.ZeroGrad);


[fmx, fpx] = BoundaryEdge.matEvaluateSurfValue( qx ); [fmy, fpy] = BoundaryEdge.matEvaluateSurfValue( qy );         
fpx = obj.matImposeNonhydroRelatedBoundaryCondition(fmx, fpx, enumNonhydroBoundaryCondition.ZeroGrad, obj.EidBoundaryType);
fpy = obj.matImposeNonhydroRelatedBoundaryCondition(fmy, fpy, enumNonhydroBoundaryCondition.ZeroGrad, obj.EidBoundaryType);
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