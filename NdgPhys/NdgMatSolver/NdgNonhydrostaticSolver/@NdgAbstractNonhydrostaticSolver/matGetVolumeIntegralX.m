function VolumeIntegralX = matGetVolumeIntegralX(obj, mesh, Variable)
%> @brief Function to calculate the volume integral in x direction
%> @details Function to calculate the volume integral in x direction
%> @param[in] mesh the mesh object
%> @param[in] Variable variable used to calculate the volume integral
%> @param[out] VolumeIntegralX the calculated volume integral in x direction

VolumeIntegralX = mesh.rx .* ( mesh.cell.Dr * Variable )...
    + mesh.sx .* (mesh.cell.Ds * Variable);

end 