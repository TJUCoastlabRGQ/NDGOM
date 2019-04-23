function VolumeIntegralY = matGetVolumeIntegralY(obj, mesh, Variable)
%> @brief Function to calculate the volume integral in y direction
%> @details Function to calculate the volume integral in y direction
%> @param[in] mesh the mesh object
%> @param[in] Variable variable used to calculate the volume integral
%> @param[out] VolumeIntegralX the calculated volume integral in y direction

VolumeIntegralY = mesh.ry .* ( mesh.cell.Dr * Variable )...
    + mesh.sy .* (mesh.cell.Ds * Variable);

end 