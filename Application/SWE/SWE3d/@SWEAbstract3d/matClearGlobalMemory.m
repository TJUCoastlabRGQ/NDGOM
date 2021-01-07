function matClearGlobalMemory( obj )
%> @brief Function to clear the global memory allocated at the begining of the simulation
%> @details
%> Function to clear the global memory allocated at the begining of the simulation
obj.advectionSolver.matClearGlobalMemory();

obj.HorizontalEddyViscositySolver.matClearGlobalMemory();

obj.PCESolver2d.matClearGlobalMemory();

obj.VerticalVelocitySolver.matClearGlobalMemory();

obj.VerticalEddyViscositySolver.matClearGlobalMemory();

end