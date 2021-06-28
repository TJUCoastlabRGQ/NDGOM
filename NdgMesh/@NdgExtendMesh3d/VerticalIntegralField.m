function [ fint3d ] = VerticalIntegralField( obj, field3d )
%VERTICALINTEGRALFIELD Summary of this function goes here
%   Detailed explanation goes here

% Np2 = obj.mesh2d.cell.Np;
% topNodeId = (obj.cell.Np - Np2 + 1):obj.cell.Np;
% Tempfint3d = zeros( obj.cell.Np, obj.K );
% fmod = obj.cell.V \ (obj.Jz .* field3d);
% 
% sk = ( (1:obj.mesh2d.K) - 1 ) * obj.Nz + obj.Nz;
% % bottom layer
% Tempfint3d(:, sk) = obj.cell.Vint * fmod(:, sk);
% cellTopField = repmat( Tempfint3d(topNodeId, sk), obj.cell.Npz, 1 );
% for n = 2:obj.Nz
%     sk = sk - 1;
%     % upper layer   bottom layer
%     Tempfint3d(:, sk) = cellTopField + obj.cell.Vint * fmod( :, sk );
%     cellTopField = repmat( Tempfint3d(topNodeId, sk), obj.cell.Npz, 1 );
% end

warning('off');
fint3d = VerticalIntegralFromBottom(struct(obj), struct(obj.cell), struct(obj.mesh2d), field3d);
warning('on');

% disp("The maximum difference is:")
% disp(max(max(fint3d - Tempfint3d)));
% disp("The minimum difference is:")
% disp(min(min(fint3d - Tempfint3d)));
end

