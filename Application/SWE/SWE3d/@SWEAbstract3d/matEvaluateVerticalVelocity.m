function [ W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d )
%MATEVALUATEVERTICALVELOCITY Summary of this function goes here
%   Detailed explanation goes here

% fphys3d(:, :, 6) = mesh3d.Extend2dField( fphys2d(:, :, 1) - fphys2d(:, :, 5) );

% HUt = mesh3d.Extend2dField( fphys2d(:, :, 2) ) - fphys3d(:, :, 5) .* fphys3d(:, :, 1);
% HVt = mesh3d.Extend2dField( fphys2d(:, :, 3) ) - fphys3d(:, :, 5) .* fphys3d(:, :, 2);

% dHUdx = mesh3d.rx .* ( mesh3d.cell.Dr * fphys3d(:, :, 1) ) + mesh3d.sx .* ( mesh3d.cell.Ds * fphys3d(:, :, 1) ) + mesh3d.tx .* ( mesh3d.cell.Dt * fphys3d(:, :, 1) );
% dHVdy = mesh3d.ry .* ( mesh3d.cell.Dr * fphys3d(:, :, 2) ) + mesh3d.sy .* ( mesh3d.cell.Ds * fphys3d(:, :, 2) ) + mesh3d.ty .* ( mesh3d.cell.Dt * fphys3d(:, :, 2) );
data = cell(1);
data{1} = fphys3d;
[ dHUHV ] = obj.matEvaluateHorizontalPartialDerivativeTerm( mesh3d,data );
% dHVdy = obj.matEvaluatePartialDerivativeTermX( mesh, fphys3d);
% fphys2d = obj.matEvaluate2dHorizonMomentum( mesh3d, ...
%         fphys2d, fphys3d );
% fphys2d(:, :, 2) = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 1) );
% fphys2d(:, :, 3) = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 2) );

% HU2D = mesh3d.Extend2dField( fphys2d(:, :, 2) );
% HV2D = mesh3d.Extend2dField( fphys2d(:, :, 3) );

data{1}(:,:,1) = mesh3d.Extend2dField( fphys2d(:, :, 2) );
data{1}(:,:,2) = mesh3d.Extend2dField( fphys2d(:, :, 3) );
[ dHU2DHV2D ] = obj.matEvaluateHorizontalPartialDerivativeTerm(mesh3d, data); 

% dHU2Ddx = mesh3d.rx .* ( mesh3d.cell.Dr * HU2D ) + mesh3d.sx .* ( mesh3d.cell.Ds * HU2D ) + mesh3d.tx .* ( mesh3d.cell.Dt * HU2D );
% dHV2Ddy = mesh3d.ry .* ( mesh3d.cell.Dr * HV2D ) + mesh3d.sy .* ( mesh3d.cell.Ds * HV2D ) + mesh3d.ty .* ( mesh3d.cell.Dt * HV2D );

W =  mesh3d.VerticalIntegralField( - dHUHV )+dHU2DHV2D.*(1+mesh3d.z) ;
end

