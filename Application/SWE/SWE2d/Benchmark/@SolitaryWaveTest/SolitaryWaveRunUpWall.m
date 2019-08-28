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








% syms x t;
% A = obj.A;
% h = obj.Depth;
% g = obj.gra;
% % x0 = 5.95;  % large amplitude
% x0 = obj.x0;
% c = sqrt(g*(A+h));
% k = sqrt(3*A/4/h^3);
% % eta =@(x,t) A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t));
% eta = A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t));
% % obj.Eta0 = eta(mesh.x,0);
% obj.Eta0 = double(subs(eta,{x,t},{mesh.x,0}));
% % U = @(x,t) c*A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t))./(A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t))+h);
% U = c*A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t))./(A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t))+h);
% difU = diff(U,x);
% % W0 = difU(mesh.x,0);
% obj.U0 = double(subs(U,{x,t},{mesh.x,0}));
% % obj.W0 = -(obj.Eta0+h).*difU(mesh.x,zeros(size(mesh.x)));
% % obj.W25 = -(obj.Eta0+h).*difU(mesh.x,25.*ones(size(mesh.x)));
% % obj.W50 = -(obj.Eta0+h).*difU(mesh.x,50.*ones(size(mesh.x)));
% W0 = zeros(size(mesh.x));
% for i = 1:size(mesh.x,1)
%     for j = 1:size(mesh.x,2)
%         W0(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),0}))*(obj.Eta0(i,j)+h);
%     end
% end
% obj.W0 = -W0;
% %% 符号变量反而计算的慢
% % obj.W0 = -double(subs(difU,{x,t},{mesh.x,0})).*(obj.Eta0+h);
% % obj.W25 = -double(subs(difU,{x,t},{mesh.x,25})).*(obj.Eta25+h);
% % obj.W50 = -double(subs(difU,{x,t},{mesh.x,50})).*(obj.Eta50+h);
end