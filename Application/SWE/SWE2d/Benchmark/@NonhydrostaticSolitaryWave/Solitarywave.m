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
        obj.H0(i,j) = double(subs(h,{x,t},{mesh.x(i,j),0}));
        obj.H5(i,j) = double(subs(h,{x,t},{mesh.x(i,j),5}));
        %         obj.Eta0(i,j) = double(subs(eta,{x,t},{mesh.x(i,j),0}));
        %         obj.Eta25(i,j) = double(subs(eta,{x,t},{mesh.x(i,j),25}));
        %         obj.Eta50(i,j) = double(subs(eta,{x,t},{mesh.x(i,j),50}));
    end
end

% % Methods from Xinhua Lu
% U = c*A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t))./(A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t))+h);
% % difU = diff(U,x);
%
% %Methods from Stelling and Zljlema
% U = c*( eta - h )/eta;
% difU = diff(U,x);


% obj.U0 = zeros(size(mesh.x));
% obj.U25 = zeros(size(mesh.x));
% obj.U50 = zeros(size(mesh.x));
% % W0 = difU(mesh.x,0);
% % obj.U0 = double(subs(U,{x,t},{mesh.x,0}));
% % obj.U25 = double(subs(U,{x,t},{mesh.x,25}));
% % obj.U50 = double(subs(U,{x,t},{mesh.x,50}));
for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.U0(i,j) = double(subs(U,{x,t},{mesh.x(i,j),0}));
        obj.U5(i,j) = double(subs(U,{x,t},{mesh.x(i,j),5}));
%         obj.U6(i,j) = double(subs(U,{x,t},{mesh.x(i,j),6}));
        %         obj.U0(i,j) = double(subs(U,{x,t},{mesh.x(i,j),0}));
        %         obj.U25(i,j) = double(subs(U,{x,t},{mesh.x(i,j),25}));
        %         obj.U50(i,j) = double(subs(U,{x,t},{mesh.x(i,j),50}));
    end
end


% % obj.W0 = -(obj.Eta0+h).*difU(mesh.x,zeros(size(mesh.x)));
% % obj.W25 = -(obj.Eta0+h).*difU(mesh.x,25.*ones(size(mesh.x)));
% % obj.W50 = -(obj.Eta0+h).*difU(mesh.x,50.*ones(size(mesh.x)));
% obj.W0 = zeros(size(mesh.x));
% W10 = zeros(size(mesh.x));
% W20 = zeros(size(mesh.x));
% obj.W25 = zeros(size(mesh.x));
% obj.W50 = zeros(size(mesh.x));
for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        %Methods from Stelling and Zljlema
        %         obj.W0(i,j) = -1 * double(subs(difU,{x,t},{mesh.x(i,j),0}))*(obj.Eta0(i,j));
        %         obj.W25(i,j) = -1 * double(subs(difU,{x,t},{mesh.x(i,j),25}))*(obj.Eta25(i,j));
        %         obj.W50(i,j) = -1 * double(subs(difU,{x,t},{mesh.x(i,j),50}))*(obj.Eta50(i,j));
        
        %         W0(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),0}))*(obj.Eta0(i,j)+h);
        %         W10(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),10}))*(obj.Eta10(i,j)+h);
        
        obj.W0(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),0} ));
        obj.W5(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),5} ));
%         obj.W6(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),6} ));
    end
end
% obj.W0 = -W0;
% obj.W25 = -W25;
% obj.W50 = -W50;
%% 符号变量反而计算的慢
% % obj.W0 = -double(subs(difU,{x,t},{mesh.x,0})).*(obj.Eta0+h);
% % obj.W25 = -double(subs(difU,{x,t},{mesh.x,25})).*(obj.Eta25+h);
% % obj.W50 = -double(subs(difU,{x,t},{mesh.x,50})).*(obj.Eta50+h);
end