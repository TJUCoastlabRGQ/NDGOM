function [ field3d ] = matVerticalRepmatFacialValue( obj, field2d )
%VERTICALINTEGRALFIELD Summary of this function goes here
%   Detailed explanation goes here

% Tempfield3d = zeros(obj.Nfp, obj.Ne);
% for i = 1:size(field2d, 2)
%     if obj.FToF(1, (i-1)*obj.Nz + 1) >= 3
%         Tempfield3d(:,(i-1)*obj.Nz + 1:i*obj.Nz) = repmat(repmat(flip(field2d(:,i)),obj.mesh.cell.Nz + 1, 1),1,obj.Nz);
%     else
%         Tempfield3d(:,(i-1)*obj.Nz + 1:i*obj.Nz) = repmat(repmat(field2d(:,i),obj.mesh.cell.Nz + 1, 1),1,obj.Nz);
%     end
% end


warning('off');
field3d = mxVerticalRepmatFacialValue(struct(obj), struct(obj.mesh.cell), field2d);
warning('on');

% disp("The maximum difference is:")
% disp(max(max(field3d - Tempfield3d)));
% disp("The minimum difference is:")
% disp(min(min(field3d - Tempfield3d)));
end