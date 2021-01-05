function matClearGlobalMemory( obj )
%> @brief Function to calculate the explicit part of the adopted time stepping method
%> @details
%> Function to calculate the explicit part of the adopted time stepping method
%> @param[in] fphys2d The two-dimensional physical field
%> @param[in] fphys The three-dimensional physical field
%> @param[in] Stage The stage of the exlicit part of the IMEXRK time stepping method 
%> @param[in] RKIndex The Runge-Kutta intermediate time step index
obj.advectionSolver.matClearGlobalMemory();

obj.HorizontalEddyViscositySolver.matClearGlobalMemory();

obj.PCESolver2d.matClearGlobalMemory();

obj.VertSolver.matClearGlobalMemory();

end