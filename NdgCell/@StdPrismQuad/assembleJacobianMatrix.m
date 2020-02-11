function [ rx, ry, rz, sx, sy, sz, tx, ty, tz, J, Jz ] = assembleJacobianMatrix( obj, x, y, z )

xr = obj.Dr * x; xs = obj.Ds * x; xt = obj.Dt * x;
yr = obj.Dr * y; ys = obj.Ds * y; yt = obj.Dt * y;
zr = obj.Dr * z; zs = obj.Ds * z; zt = obj.Dt * z;
J = xr .* (ys .* zt - zs .* yt ) - yr .* ( xs .* zt - zs .* xt ) ...
    + zr .* (xs .* yt - ys .* xt);

rx = ( ys .* zt - zs .* yt ) ./ J;
sx = - (yr .* zt - zr .* yt) ./ J;
tx = ( yr .* zs - zr .* ys ) ./ J;
ry = - ( xs .* zt - zs .* xt ) ./ J;
sy = ( xr .* zt - zr .* xt ) ./ J;
ty = - ( xr .* zs - zr .* xs ) ./ J;
rz = ( xs .* yt - ys .* xt ) ./ J;
sz = - ( xr .* yt - yr .* xt ) ./ J;
tz = ( xr .* ys - yr .* xs ) ./ J;

Jz = zt;
end