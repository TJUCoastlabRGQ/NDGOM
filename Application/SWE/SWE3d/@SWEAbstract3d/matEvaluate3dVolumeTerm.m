function Volume_frhs3d = matEvaluate3dVolumeTerm( obj, mesh3d, fphys3d )
%> @brief Function used to calculate the volume integration part of the model
%> 
%> More detailed description.
%>
%> @param mesh3d The three dimensional mesh object
%> @param fphys3d The three dimensional physical field
%>
%> @retval Volume_frhs3d The three dimensional volume integration part


HUX = fphys3d(:, :, 1).^2./fphys3d(:, :, 6) + 1/2 * obj.gra * fphys3d(:, :, 6).^2;   % Hu^2 + 1/2gH^2 
HUY = fphys3d(:, :, 1) .* fphys3d(:, :, 2)./fphys3d(:, :, 6);   % Huv

HUZ = fphys3d(:, :, 1) .* fphys3d(:, :, 3)./fphys3d(:, :, 6) - fphys3d(:, :, 4);  %u*omega - taux

HVX = fphys3d(:, :, 1) .* fphys3d(:, :, 2)./fphys3d(:, :, 6);   % Huv
HVY = fphys3d(:, :, 2).^2./fphys3d(:, :, 6) + 1/2 * obj.gra * fphys3d(:, :, 6).^2;   % Hv^2 + 1/2gH^2 
HVZ = fphys3d(:, :, 2) .* fphys3d(:, :, 3)./fphys3d(:, :, 6) - fphys3d(:, :, 5);  %v*omega - tauy


Volume_frhs3d(:, :, 1) =...
    - ( mesh3d.rx .* (mesh3d.cell.Dr * HUX) + mesh3d.sx .* (mesh3d.cell.Ds * HUX) + mesh3d.tx .* (mesh3d.cell.Dt * HUX) + ...
        mesh3d.ry .* (mesh3d.cell.Dr * HUY) + mesh3d.sy .* (mesh3d.cell.Ds * HUY) + mesh3d.ty .* (mesh3d.cell.Dt * HUY) +...
        mesh3d.rz .* (mesh3d.cell.Dr * HUZ) + mesh3d.sz .* (mesh3d.cell.Ds * HUZ) + mesh3d.tz .* (mesh3d.cell.Dt * HUZ) ) ;

Volume_frhs3d(:, :, 2) =...
    - ( mesh3d.rx .* (mesh3d.cell.Dr * HVX) + mesh3d.sx .* (mesh3d.cell.Ds * HVX) + mesh3d.tx .* (mesh3d.cell.Dt * HVX) +...
        mesh3d.ry .* (mesh3d.cell.Dr * HVY) + mesh3d.sy .* (mesh3d.cell.Ds * HVY) + mesh3d.ty .* (mesh3d.cell.Dt * HVY) + ...
        mesh3d.rz .* (mesh3d.cell.Dr * HVZ) + mesh3d.sz .* (mesh3d.cell.Ds * HVZ) + mesh3d.tz .* (mesh3d.cell.Dt * HVZ) );

end