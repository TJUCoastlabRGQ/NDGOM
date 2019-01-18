function [termx, termy] = matCalculateLDGAuxialaryTerm( obj, mesh, BoundaryEdge, InnerEdge, variable )
%> @brief Function to calculate variable with the upwinded flux to weakly impose the continuity condition
%> @details
%> Function to calculate variable with the upwinded flux to weakly impose the continuity condition
%> @param[in] mesh The mesh object
%> @param[in] BoundaryEdge The boundaryEdge object
%> @param[in] InnerEdge The innerEdge object
%> @param[in] Variable The variable used to calculate the returned value
%> @param[out] termx the calculated characteristic matrix in x direction
%> @param[out] termy the calculated characteristic matrix in y direction
[fm, fp] = InnerEdge.matEvaluateSurfValue( variable );       
[fm, fp] = obj.matGetFaceValue(fm, fp, enumNonhydroBoundaryCondition.Zero);
%< Inner edge contribution
fluxMX = InnerEdge.nx.*fm; fluxMY = InnerEdge.ny.*fm;
fluxPX = InnerEdge.nx.*fp; fluxPY = InnerEdge.ny.*fp; 

Beta = [1, 1];
[penaltyX, penaltyY] = obj.matEvaluateInnerEdgeScalarJumpTerm( InnerEdge, variable, enumNonhydroBoundaryCondition.Zero);

fluxSX = InnerEdge.nx .* ((fm + fp)./2 - Beta(1) * penaltyX - Beta(2) * penaltyY);
fluxSY = InnerEdge.ny .* ((fm + fp)./2 - Beta(1) * penaltyX - Beta(2) * penaltyY);

termx = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMX, fluxPX, fluxSX );
termy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMY, fluxPY, fluxSY );

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( variable );        
fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, enumNonhydroBoundaryCondition.Zero, obj.EidBoundaryType);
%< Boundary edge contribution
fluxMX = BoundaryEdge.nx.*fm; fluxMY = BoundaryEdge.ny.*fm;


[penaltyX, penaltyY] = obj.matEvaluateBoundaryEdgeScalarJumpTerm( BoundaryEdge, variable, enumNonhydroBoundaryCondition.Zero);
fluxSX = BoundaryEdge.nx .* ((fm + fp)./2 - Beta(1) * penaltyX - Beta(2) * penaltyY);
fluxSY = BoundaryEdge.ny .* ((fm + fp)./2 - Beta(1) * penaltyX - Beta(2) * penaltyY);

termx = - termx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
termy = - termy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);

termx = termx + mesh.rx .* (mesh.cell.Dr * cell2mat(variable))...
    + mesh.sx .* (mesh.cell.Ds * cell2mat(variable));
termy = termy + mesh.ry .* (mesh.cell.Dr * cell2mat(variable))...
    + mesh.sy .* (mesh.cell.Ds * cell2mat(variable));



end