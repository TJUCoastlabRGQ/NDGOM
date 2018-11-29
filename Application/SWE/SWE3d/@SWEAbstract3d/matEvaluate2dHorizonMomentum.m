function fphys2d = matEvaluate2dHorizonMomentum( obj, mesh3d, fphys2d, fphys3d )

% evaluate 2d averaged velocity
U = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 1) );
V = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 2) );

%fphys2d(:, :, 4) = fphys2d(:, :, 1) - fphys2d(:, :, 5);
fphys2d(:, :, 2) = fphys2d(:, :, 4) .* U;
fphys2d(:, :, 3) = fphys2d(:, :, 4) .* V;

end