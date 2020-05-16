function [ field2d ] = VerticalColumnIntegralField( obj, field3d )
%VERTICALINTEGRALFIELD Summary of this function goes here
%   Detailed explanation goes here

Np2 = size(obj.V1d,1);
K2d = size(field3d,2)/obj.Nz;

fmod = obj.V2d \ (obj.Jz .* field3d);

field2d = zeros( Np2, K2d );

sk = ( (1:K2d) - 1 ) * obj.Nz + 1;
for n = 1:obj.Nz
    field2d = field2d + fmod( 1:Np2, sk );
    sk = sk + 1;
end

% the mon coefficient need to devide the \varphi_{z,0}
% field2d = obj.mesh2d.cell.V * ( field2d ./ obj.mesh2d.cell.V(:, 1) );
% field2d = obj.mesh2d.cell.V * bsxfun(@rdivide, field2d, obj.mesh2d.cell.V(:, 1) );
field2d = obj.V1d * field2d / (sqrt(2)/2);

% the date need to be flipped, because the order of face point between 3d
% face and 2d line is not the same
Index = ( obj.FToF(1,1:obj.Nz:end) >= 3 );
field2d(:,Index) = flip(field2d(:,Index));
end
