function Solitarywave(obj, mesh)


% if count(py.sys.path,'') == 0
%     insert(py.sys.path,int32(0),'');
% end
% 
% Module = py.importlib.import_module('sympy');
% 
% x  = Module.Symbol('x');t  = Module.Symbol('t');
% 
% % f = Module.Symbol('f', cls=Function );
% f = 2*x+t;

% f.evalf( Module.subs={x:1, t:2})

syms x t;

H0 = 1;
d = H0; % Constant produce during the derivation, and is taken to be equal to the water depth
A = 0.2;
l = H0*sqrt((A+H0)/A);
C0 = l/d*sqrt(obj.gra * H0^3/(l^2-H0^2));
% Ce = sqrt(obj.gra *(A+H0));
h = H0 + A * ( sech( (x-C0*t)/l ) )^2;
U = C0 * (1 - d/h);
% W = -( A * Ce * d )/( l * h ) * sech((x - Ce * t)/l)*( diff(sech((x - Ce * t)/l),x) + diff(sech((x - Ce * t)/l),t) );
% W = -( A * Ce * d )/( l * h ) * sech((x - Ce * t)/l)*( diff(sech((x - Ce * t)/l),t) );
W = -( A * C0 * d )/( l * h ) * sech( (x - C0 * t)/l )*( diff( sech( (x - C0 * t)/l ), x ) * l  );
P = A*C0^2*d^2/2/l^2/h^2*((2*H0-h)*(diff(sech((x - C0 * t)/l),x) * l )^2+...
    h * sech((x-C0*t)/l) * diff(diff(sech((x - C0 * t)/l),x) * l, x ) * l);
% f = diff(h*U,x) - U * diff(h,x) + 2*W;



for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.P0(i,j) = double(subs(P,{x,t},{mesh.x(i,j),0}));
        obj.P5(i,j) = double(subs(P,{x,t},{mesh.x(i,j),5}));
%         obj.P8(i,j) = double(subs(P,{x,t},{mesh.x(i,j),8}));
%         obj.P12(i,j) = double(subs(P,{x,t},{mesh.x(i,j),12}));
%         obj.P16(i,j) = double(subs(P,{x,t},{mesh.x(i,j),16}));
        obj.H0(i,j) = double(subs(h,{x,t},{mesh.x(i,j),0}));
        obj.H5(i,j) = double(subs(h,{x,t},{mesh.x(i,j),5}));
%         obj.H8(i,j) = double(subs(h,{x,t},{mesh.x(i,j),8}));
%         obj.H12(i,j) = double(subs(h,{x,t},{mesh.x(i,j),12}));
%         obj.H16(i,j) = double(subs(h,{x,t},{mesh.x(i,j),16}));
    end
end

for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.U0(i,j) = double(subs(U,{x,t},{mesh.x(i,j),0}));
        obj.U5(i,j) = double(subs(U,{x,t},{mesh.x(i,j),5}));
%         obj.U8(i,j) = double(subs(U,{x,t},{mesh.x(i,j),8}));
%         obj.U12(i,j) = double(subs(U,{x,t},{mesh.x(i,j),12}));
%         obj.U16(i,j) = double(subs(U,{x,t},{mesh.x(i,j),16}));
    end
end

for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.W0(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),0} ));
        obj.W5(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),5} ));
%         obj.W8(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),8} ));
%         obj.W12(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),12} ));
%         obj.W16(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),16} ));
    end
end
end