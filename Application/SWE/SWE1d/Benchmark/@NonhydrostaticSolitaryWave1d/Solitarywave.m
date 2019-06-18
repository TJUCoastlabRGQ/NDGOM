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
        obj.U0(i,j) = double(subs(U,{x,t},{mesh.x(i,j),0}));
        obj.W0(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),0} ));
        obj.H0(i,j) = double(subs(h,{x,t},{mesh.x(i,j),0}));
    end
end
       obj.U5 = zeros(size(mesh.x));
       obj.H5 = zeros(size(mesh.x));
       obj.P5 = zeros(size(mesh.x));
       obj.W5 = zeros(size(mesh.x));
for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
%         T = ( obj.xlim(2) - obj.xlim(1) ) / C0;
%         tempt = mod(C0*5, obj.xlim(2) - obj.xlim(1) )/C0;
%         Index = mesh.x + mod(C0*5, obj.xlim(2) - obj.xlim(1) ) > obj.xlim(2);
%         obj.U5(Index) = double(subs(U,{x,t},{mesh.x(Index), - (T - tempt)}));
%         obj.H5(Index) = double(subs(h,{x,t},{mesh.x(Index), - (T - tempt)}));
%         obj.P5(Index) = double(subs(P,{x,t},{mesh.x(Index), - (T - tempt)}));
%         obj.W5(Index) = double(subs(W,{x,t},{mesh.x(Index), - (T - tempt)} ));
%         
%         Index = mesh.x + mod(C0*5, obj.xlim(2) - obj.xlim(1) ) <= obj.xlim(2);
% 
%         obj.U5(Index) = double(subs(U,{x,t},{mesh.x(Index),  tempt}));
%         obj.H5(Index) = double(subs(h,{x,t},{mesh.x(Index),  tempt}));
%         obj.P5(Index) = double(subs(P,{x,t},{mesh.x(Index),  tempt}));
%         obj.W5(Index) = double(subs(W,{x,t},{mesh.x(Index),  tempt}));
        
        obj.U5(i,j) = double(subs(U,{x,t},{mesh.x(i,j),5}));
        obj.H5(i,j) = double(subs(h,{x,t},{mesh.x(i,j),5}));
        obj.P5(i,j) = double(subs(P,{x,t},{mesh.x(i,j),5}));
        obj.W5(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),5} ));
    end
end
end