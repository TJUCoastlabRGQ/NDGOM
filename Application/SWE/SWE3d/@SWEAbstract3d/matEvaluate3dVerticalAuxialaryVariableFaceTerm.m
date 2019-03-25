function [ AuxialaryVariableFace_rhs3d ] = matEvaluate3dVerticalAuxialaryVariableFaceTerm( obj, mesh3d, fphys3d )
%> @brief Function used to calculate the surface contribution of the auxialary variable to the RHS
%> 
%> More detailed description.
%>
%> For this case, we only consider the vertical auxialary variable contribution to the
%> right hand side, and the contribution due to the vertical momentum term is
%> not considered here. For this considered part, the Nueman boundary
%> condition is considered and  '$\mathbf \tau^+ = -\mathbf \tau^- +
%> 2\mathbf \tau_s$' at '$\omega = 0$', and '$\mathbf \tau^+ = -\mathbf \tau^- +
%> 2 C_f\mathbf u_{xy}(-1)||\mathbf u_{xy}(-1)||$' at '$\omega = -1$'
%>
%> @param mesh3d The three dimensional mesh object
%> @param fphys3d The three dimensional physical field
%>
%> @retval AuxialaryVariableFace_rhs3d The right hand side due to the auxialary variable in vertical direction

edge = mesh3d.BottomEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );
fluxM(:,:,1) = fm(:,:,4); fluxM(:,:,2) = fm(:,:,5);
fluxP(:,:,1) = fp(:,:,4); fluxP(:,:,2) = fp(:,:,5);
AuxialaryVariableFace_rhs3d = edge.matEvaluateStrongFormEdgeCentralRHS( fluxM, fluxP );
clear fluxM fluxP;

edge = mesh3d.SurfaceBoundaryEdge;
[ fm, ~ ] = edge.matEvaluateSurfValue( fphys3d );
fluxM(:,:,1) = fm(:,:,4); 
fluxM(:,:,2) = fm(:,:,5);
fluxS(:,:,1) = obj.Taux; fluxS(:,:,2) = obj.Tauy;

% At surface, the unit vector is 1
AuxialaryVariableFace_rhs3d = AuxialaryVariableFace_rhs3d - ...
    edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );
clear fluxM fluxS;

edge = mesh3d.BottomBoundaryEdge;
[ fm, ~ ] = edge.matEvaluateSurfValue( fphys3d );
fluxM(:,:,1) = fm(:,:,4); fluxM(:,:,2) = fm(:,:,5);
u = fm(:,:,1)./fm(:,:,6); v = fm(:,:,2)./fm(:,:,6); 
Velocity = sqrt( u.^2 + v.^2 );
fluxS(:,:,1) = obj.Cf .* u .* Velocity;
fluxS(:,:,2) = obj.Cf .* v .* Velocity;

AuxialaryVariableFace_rhs3d = AuxialaryVariableFace_rhs3d + ...
    edge.matEvaluateStrongFormEdgeRHS( fluxM, fluxS );

end