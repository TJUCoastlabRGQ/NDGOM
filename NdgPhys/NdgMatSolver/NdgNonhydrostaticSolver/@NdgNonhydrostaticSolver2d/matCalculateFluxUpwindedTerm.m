function [termx, termy] = matCalculateFluxUpwindedTerm(obj, mesh, BoundaryEdge, InnerEdge, Variable, UpWindedFlag, DownWindedFlag, eidtype)
%> @brief Function to calculate variable with the upwinded flux to weakly impose the continuity condition
%> @details
%> Function to calculate variable with the upwinded flux to weakly impose the continuity condition
%> @param[in] mesh The mesh object
%> @param[in] BoundaryEdge The boundaryEdge object
%> @param[in] InnerEdge The innerEdge object
%> @param[in] Variable The variable used to calculate the returned value
%> @param[in] eidtype The boundary condition to be imposed at the wet-dry interface and all the clamped-type boundary 
%> @param[out] termx the calculated characteristic matrix in x direction
%> @param[out] termy the calculated characteristic matrix in y direction
[fm, fp] = InnerEdge.matEvaluateSurfValue( Variable );       
[fm, fp] = obj.matGetFaceValue(fm, fp, eidtype);
%< Inner edge contribution
fluxMX = InnerEdge.nx.*fm; fluxMY = InnerEdge.ny.*fm;
fluxPX = InnerEdge.nx.*fp; fluxPY = InnerEdge.ny.*fp; 
% fluxSx = InnerEdge.nx .* (1 + sign(InnerEdge.nx .* fm + InnerEdge.ny .* fm))./2 .* fm + InnerEdge.nx .*  (1 + sign( -InnerEdge.nx .* fm - InnerEdge.ny .* fm))./2 .* fp;
% fluxSy = InnerEdge.ny .* (1 + sign(InnerEdge.nx .* fm + InnerEdge.ny .* fm))./2 .* fm + InnerEdge.ny .*  (1 + sign( -InnerEdge.nx .* fm - InnerEdge.ny .* fm))./2 .* fp;
% fluxSx = InnerEdge.nx .* fm;
% fluxSy = InnerEdge.ny .* fm;
% fluxSx = InnerEdge.nx .* (1 + sign(InnerEdge.nx  + InnerEdge.ny ))./2 .* fm + InnerEdge.nx .*  (1 + sign( -InnerEdge.nx  - InnerEdge.ny ))./2 .* fp;
% fluxSy = InnerEdge.ny .* (1 + sign(InnerEdge.nx  + InnerEdge.ny ))./2 .* fm + InnerEdge.ny .*  (1 + sign( -InnerEdge.nx  - InnerEdge.ny ))./2 .* fp;
fluxSx = UpWindedFlag .* InnerEdge.nx .* fm - DownWindedFlag .* InnerEdge.nx .* fp;
fluxSy = UpWindedFlag .* InnerEdge.ny .* fm - DownWindedFlag .* InnerEdge.ny .* fp;

termx = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMX, fluxPX, fluxSx );
termy = InnerEdge.matEvaluateStrongFromEdgeRHS( fluxMY, fluxPY, fluxSy );

[fm, fp] = BoundaryEdge.matEvaluateSurfValue( Variable );        
fp = obj.matImposeNonhydroRelatedBoundaryCondition(fm, fp, eidtype, obj.EidBoundaryType);
%< Boundary edge contribution
fluxMX = BoundaryEdge.nx.*fm; fluxMY = BoundaryEdge.ny.*fm;
% fluxSX = BoundaryEdge.nx .* (1 + sign(BoundaryEdge.nx .* fm + BoundaryEdge.ny .* fm))./2 .* fm + BoundaryEdge.nx .*  (1 + sign( -BoundaryEdge.nx .* fm - BoundaryEdge.ny .* fm))./2 .* fp;
% fluxSY = BoundaryEdge.ny .* (1 + sign(BoundaryEdge.nx .* fm + BoundaryEdge.ny .* fm))./2 .* fm + BoundaryEdge.ny .*  (1 + sign( -BoundaryEdge.nx .* fm - BoundaryEdge.ny .* fm))./2 .* fp;
% fluxSX = BoundaryEdge.nx .* fm;
% fluxSY = BoundaryEdge.ny .* fm;
% fluxSX = BoundaryEdge.nx .* (1 + sign(BoundaryEdge.nx + BoundaryEdge.ny))./2 .* fm + BoundaryEdge.nx .*  (1 + sign( -BoundaryEdge.nx - BoundaryEdge.ny))./2 .* fp;
% fluxSY = BoundaryEdge.ny .* (1 + sign(BoundaryEdge.nx + BoundaryEdge.ny))./2 .* fm + BoundaryEdge.ny .*  (1 + sign( -BoundaryEdge.nx - BoundaryEdge.ny))./2 .* fp;
fluxSX = BoundaryEdge.nx .* fm;
fluxSY = BoundaryEdge.ny .* fm;

termx = - termx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX, fluxSX);
termy = - termy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY, fluxSY);

termx = termx + mesh.rx .* (mesh.cell.Dr * cell2mat(Variable))...
    + mesh.sx .* (mesh.cell.Ds * cell2mat(Variable));
termy = termy + mesh.ry .* (mesh.cell.Dr * cell2mat(Variable))...
    + mesh.sy .* (mesh.cell.Ds * cell2mat(Variable));

end