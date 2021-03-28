function matCalculateExplicitRHSTerm( obj, fphys, time, Stage, RKindex )
%> @brief Function to calculate the explicit part of the adopted time stepping method
%> @details
%> Function to calculate the explicit part of the adopted time stepping method
%> @param[in] fphys The three-dimensional physical field
%> @param[in] time The local time of the current calculation
obj.advectionSolver.evaluateAdvectionRHS( obj, fphys );
obj.HorizontalEddyViscositySolver.matEvaluateDiffRHS(obj, fphys);
obj.matEvaluateSourceTerm( time );
obj.ExplicitRHS(:,:,RKindex:Stage:end) = obj.frhs{1}(:,:,:);
end