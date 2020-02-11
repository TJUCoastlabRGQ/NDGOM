function SolitaryWaveRunUpWall(obj, mesh)


syms x t;

H0 = 0.218;
d = H0; % Constant produce during the derivation, and is taken to be equal to the water depth
A = obj.A;
x0 = obj.x0;
l = H0*sqrt((A+H0)/A);
C0 = l/d*sqrt(obj.gra * H0^3/(l^2-H0^2));
% Ce = sqrt(obj.gra *(A+H0));
h = H0 + A * ( sech( (x - x0 - C0*t)/l ) )^2;
U = C0 * (1 - d/h);
% W = -( A * Ce * d )/( l * h ) * sech((x - Ce * t)/l)*( diff(sech((x - Ce * t)/l),x) + diff(sech((x - Ce * t)/l),t) );
% W = -( A * Ce * d )/( l * h ) * sech((x - Ce * t)/l)*( diff(sech((x - Ce * t)/l),t) );
W = -( A * C0 * d )/( l * h ) * sech( (x - x0 - C0 * t)/l )*( diff( sech( (x  - x0 - C0 * t)/l ), x ) * l  );
% f = diff(h*U,x) - U * diff(h,x) + 2*W;



for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.Eta0(i,j) = double(subs(h,{x,t},{mesh.x(i,j),0}));
    end
end

for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.U0(i,j) = double(subs(U,{x,t},{mesh.x(i,j),0}));
    end
end

for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.W0(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),0} ));
    end
end

end