function Solitarywave(obj, mesh)
syms x t;

H0 = 1;
d = H0; % Constant produce during the derivation, and is taken to be equal to the water depth
A = 0.2;
l = H0*sqrt((A+H0)/A);
Ce = l/d*sqrt(obj.gra * H0^3/(l^2-H0^2));

H = H0 + A * ( sech((x-Ce*t)/l) )^2;
U = Ce * (1 - d/H);
W = -( A * Ce * d )/( l * H ) * sech((x - Ce * t)/l)*diff(sech((x - Ce * t)/l),x);
P = A*Ce^2*d^2/2/l^2/H^2*((2*H0-H)*(diff(sech((x - Ce * t)/l),x))^2+...
    H * sech((x-Ce*t)/l) * diff(diff(sech((x - Ce * t)/l),x), x ));
% A = obj.A;
% h = obj.Depth;
% g = obj.gra;
% x0 = 200;
% c = sqrt(g*(A+h));
% k = sqrt(3*A/4/h^3);
% Methods from Xinhua Lu
% eta = A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t));

%Methods from Stelling and Zljlema
% eta = h + 4*A*(exp(-sqrt(3*A/h/h/(h+A))*(x-c*t - x0)))/(1 + exp(-sqrt(3*A/h/h/(h+A))*(x-c*t - x0)))^2;

% eta = 4*A*exp(-sqrt(3*A/(h^2*(h+A)))*(x-x0-c*t))./(1+exp(-sqrt(3*A/(h^2*(h+A)))*(x-x0-c*t)))^2;
% obj.Eta0 = eta(mesh.x,0);
% obj.Eta0 = double(subs(eta,{x,t},{mesh.x,0}));
% obj.Eta25 = double(subs(eta,{x,t},{mesh.x,25}));
% obj.Eta50 = double(subs(eta,{x,t},{mesh.x,50}));

% obj.Eta0 = zeros(size(mesh.x));
% obj.Eta25 = zeros(size(mesh.x));
% obj.Eta50 = zeros(size(mesh.x));

for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.P0(i,j) = double(subs(P,{x,t},{mesh.x(i,j),0}));
        obj.P6(i,j) = double(subs(P,{x,t},{mesh.x(i,j),6}));
        obj.H0(i,j) = double(subs(H,{x,t},{mesh.x(i,j),0}));
        obj.H3(i,j) = double(subs(H,{x,t},{mesh.x(i,j),3}));
        obj.H6(i,j) = double(subs(H,{x,t},{mesh.x(i,j),6}));
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
        obj.U3(i,j) = double(subs(U,{x,t},{mesh.x(i,j),3}));
        obj.U6(i,j) = double(subs(U,{x,t},{mesh.x(i,j),6}));
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
        obj.W3(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),3} ));
        obj.W6(i,j) = double(subs(W, {x,t}, {mesh.x(i,j),6} ));
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