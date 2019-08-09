function Solitarywave(obj, mesh)

syms x t;

H0 = 0.32;
d = H0;
A = H0 * obj.Ratio;
l = H0*sqrt((A+H0)/A);
C0 = l/d*sqrt(obj.gra * H0^3/(l^2-H0^2));
h = H0 + A * ( sech( (x-C0*t)/l ) )^2;
U = C0 * (1 - d/h);
W = -( A * C0 * d )/( l * h ) * sech( (x - C0 * t)/l )*( diff( sech( (x - C0 * t)/l ), x )  * l  );
% P = A*C0^2*d^2/2/l^2/h^2*((2*H0-h)*(diff(sech((x - C0 * t)/l),x) * l )^2+...
%     h * sech((x-C0*t)/l) * diff(diff(sech((x - C0 * t)/l),x) * l, x ) * l);



for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
%         obj.P0(i,j) = double(subs(P,{x,t},{mesh.x(i,j),0}));
        obj.H0(i,j) = double(subs(h,{x,t},{mesh.x(i,j),0}));
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