function  [ VolumeIntegralX ] = matVolumeIntegral( obj, mesh, VariableX )
%> @brief Function to calculate the volume integral in y direction
%> @details Function to calculate the volume integral in y direction
%> @param[in] mesh the mesh object
%> @param[in] Variable variable used to calculate the volume integral
%> @param[out] VolumeIntegralX the calculated volume integral in y direction

VolumeIntegralX = mesh.rx .* ( mesh.cell.Dr * VariableX );

end