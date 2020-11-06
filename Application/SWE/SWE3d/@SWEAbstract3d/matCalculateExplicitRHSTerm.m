function matCalculateExplicitRHSTerm( obj, fphys2d, fphys, Stage, RKIndex)
%> @brief Function to calculate the explicit part of the adopted time stepping method
%> @details
%> Function to calculate the explicit part of the adopted time stepping method
%> @param[in] fphys2d The two-dimensional physical field
%> @param[in] fphys The three-dimensional physical field
%> @param[in] Stage The stage of the exlicit part of the IMEXRK time stepping method 
%> @param[in] RKIndex The Runge-Kutta intermediate time step index
obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
% obj.HorizontalEddyViscositySolver.matEvaluateDiffRHS(obj, fphys{1});
obj.PCESolver2d.evaluateAdvectionRHS(obj, fphys2d );
obj.matEvaluateSourceTerm( fphys );
obj.ExplicitRHS2d(:,:,RKIndex) = obj.frhs2d{1}(:,:,1);
obj.ExplicitRHS(:,:,RKIndex:Stage:end) = obj.frhs{1}(:,:,:);

end