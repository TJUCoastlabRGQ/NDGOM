function frhs2d_VolumeTerm = matEvaluate2dHorizonPCEVolumeTerm( obj, mesh2d, fphys2d )
%> @brief Function used to calculate the two dimentional PCE volume term
%> 
%> More detailed description.
%>
%> @param mesh2d The two dimensional mesh object
%> @param fphys2d The two dimensional physical field
%>
%> @fphys2d frhs2d_VolumeTerm Contribution to the RHS of the two dimensional PCE from the volume integration term

frhs2d_VolumeTerm(:, :, 1) = -( ...
    mesh2d.rx .* ( mesh2d.cell.Dr * fphys2d(:, :, 2) ) + ...
    mesh2d.sx .* ( mesh2d.cell.Ds * fphys2d(:, :, 2) ) + ...
    mesh2d.ry .* ( mesh2d.cell.Dr * fphys2d(:, :, 3) ) + ...
    mesh2d.sy .* ( mesh2d.cell.Ds * fphys2d(:, :, 3) ) );

end