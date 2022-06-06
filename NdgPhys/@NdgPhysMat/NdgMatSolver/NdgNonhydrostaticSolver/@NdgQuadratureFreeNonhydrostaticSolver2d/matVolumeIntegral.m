function  [ VolumeIntegralX , VolumeIntegralY ] = matVolumeIntegral(obj, mesh, VariableX, VariableY)
%> @brief Function to calculate the volume integral in y direction
%> @details Function to calculate the volume integral in y direction
%> @param[in] mesh the mesh object
%> @param[in] Variable variable used to calculate the volume integral
%> @param[out] VolumeIntegralX the calculated volume integral in y direction

VolumeIntegralX = mesh.rx .* ( mesh.cell.Dr * VariableX )...
    + mesh.sx .* (mesh.cell.Ds * VariableX);

VolumeIntegralY = mesh.ry .* ( mesh.cell.Dr * VariableY )...
    + mesh.sy .* (mesh.cell.Ds * VariableY);

end