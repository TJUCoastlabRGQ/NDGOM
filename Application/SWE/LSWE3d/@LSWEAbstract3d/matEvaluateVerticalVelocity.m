function [ W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d )
%MATEVALUATEVERTICALVELOCITY Summary of this function goes here
%   Detailed explanation goes here

% % fphys3d(:, :, 6) = mesh3d.Extend2dField( fphys2d(:, :, 1) - fphys2d(:, :, 5) );
% 
% % HUt = mesh3d.Extend2dField( fphys2d(:, :, 2) ) - fphys3d(:, :, 5) .* fphys3d(:, :, 1);
% % HVt = mesh3d.Extend2dField( fphys2d(:, :, 3) ) - fphys3d(:, :, 5) .* fphys3d(:, :, 2);
% 
% HUt = - fphys3d(:, :, 6) .* fphys3d(:, :, 1);
% HVt = - fphys3d(:, :, 6) .* fphys3d(:, :, 2);
% 
% dHUdx = mesh3d.rx .* ( mesh3d.cell.Dr * HUt ) + mesh3d.sx .* ( mesh3d.cell.Ds * HUt );
% dHVdy = mesh3d.ry .* ( mesh3d.cell.Dr * HVt ) + mesh3d.sy .* ( mesh3d.cell.Ds * HVt );
% 
% HU2D = mesh3d.Extend2dField( fphys2d(:, :, 2) );
% HV2D = mesh3d.Extend2dField( fphys2d(:, :, 3) );
% dHU2Ddx = mesh3d.rx .* ( mesh3d.cell.Dr * HU2D ) + mesh3d.sx .* ( mesh3d.cell.Ds * HU2D );
% dHV2Ddy = mesh3d.ry .* ( mesh3d.cell.Dr * HV2D ) + mesh3d.sy .* ( mesh3d.cell.Ds * HV2D );
% 
% % W = mesh3d.VerticalIntegralField( dHUdx + dHVdy );
% W = mesh3d.VerticalIntegralField( dHUdx + dHVdy )+dHU2Ddx.*(1+mesh3d.z)+dHV2Ddy.*(1+mesh3d.z);

data = cell(1);
data{1} = fphys3d;
data{1}(:,:,1) = mesh3d.Extend2dField( fphys2d(:, :, 4) ) .* data{1}(:,:,1);
data{1}(:,:,2) = mesh3d.Extend2dField( fphys2d(:, :, 4) ) .* data{1}(:,:,2);
[ dHUdx, dHVdy ] = obj.matEvaluatePartialDerivativeTermX( mesh3d,data );
% dHVdy = obj.matEvaluatePartialDerivativeTermX( mesh, fphys3d);
% fphys2d = obj.matEvaluate2dHorizonMomentum( mesh3d, ...
%         fphys2d, fphys3d );

HU2D = mesh3d.Extend2dField( fphys2d(:, :, 2) );
HV2D = mesh3d.Extend2dField( fphys2d(:, :, 3) );
dHU2Ddx = mesh3d.rx .* ( mesh3d.cell.Dr * HU2D ) + mesh3d.sx .* ( mesh3d.cell.Ds * HU2D ) + mesh3d.tx .* ( mesh3d.cell.Dt * HU2D );
dHV2Ddy = mesh3d.ry .* ( mesh3d.cell.Dr * HV2D ) + mesh3d.sy .* ( mesh3d.cell.Ds * HV2D ) + mesh3d.ty .* ( mesh3d.cell.Dt * HV2D );

W = mesh3d.VerticalIntegralField( - dHUdx - dHVdy )+dHU2Ddx.*(1+mesh3d.z)+dHV2Ddy.*(1+mesh3d.z);
end

