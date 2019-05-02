function Solitarywave(obj, mesh)
syms x t;
A = obj.A;
h = obj.Depth;
g = obj.gra;
x0 = 200;
c = sqrt(g*(A+h));
k = sqrt(3*A/4/h^3);
% Methods from Xinhua Lu
eta = A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t));

%Methods from Stelling and Zljlema
% eta = h + 4*A*(exp(-sqrt(3*A/h/h/(h+A))*(x-c*t - x0)))/(1 + exp(-sqrt(3*A/h/h/(h+A))*(x-c*t - x0)))^2;

% eta = 4*A*exp(-sqrt(3*A/(h^2*(h+A)))*(x-x0-c*t))./(1+exp(-sqrt(3*A/(h^2*(h+A)))*(x-x0-c*t)))^2;
% obj.Eta0 = eta(mesh.x,0);
% obj.Eta0 = double(subs(eta,{x,t},{mesh.x,0}));
% obj.Eta25 = double(subs(eta,{x,t},{mesh.x,25}));
% obj.Eta50 = double(subs(eta,{x,t},{mesh.x,50}));

for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.Eta0(i,j) = double(subs(eta,{x,t},{mesh.x(i,j),0}));
        obj.Eta1(i,j) = double(subs(eta,{x,t},{mesh.x(i,j),1}));
        obj.Eta2(i,j) = double(subs(eta,{x,t},{mesh.x(i,j),2}));
%         obj.Eta25(i,j) = double(subs(eta,{x,t},{mesh.x(i,j),25}));
%         obj.Eta50(i,j) = double(subs(eta,{x,t},{mesh.x(i,j),50}));
    end
end

% Methods from Xinhua Lu
U = c*A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t))./(A./cosh(k*(x-x0-c*t))./cosh(k*(x-x0-c*t))+h);
% difU = diff(U,x);

%Methods from Stelling and Zljlema
% U = c*( eta - h )/eta;
difU = diff(U,x);


% % W0 = difU(mesh.x,0);
% % obj.U0 = double(subs(U,{x,t},{mesh.x,0}));
% % obj.U25 = double(subs(U,{x,t},{mesh.x,25}));
% % obj.U50 = double(subs(U,{x,t},{mesh.x,50}));
for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        obj.U0(i,j) = double(subs(U,{x,t},{mesh.x(i,j),0}));
        obj.U1(i,j) = double(subs(U,{x,t},{mesh.x(i,j),1}));
        obj.U2(i,j) = double(subs(U,{x,t},{mesh.x(i,j),2}));
%         obj.U25(i,j) = double(subs(U,{x,t},{mesh.x(i,j),25}));
%         obj.U50(i,j) = double(subs(U,{x,t},{mesh.x(i,j),50}));
    end
end


% % obj.W0 = -(obj.Eta0+h).*difU(mesh.x,zeros(size(mesh.x)));
% % obj.W25 = -(obj.Eta0+h).*difU(mesh.x,25.*ones(size(mesh.x)));
% % obj.W50 = -(obj.Eta0+h).*difU(mesh.x,50.*ones(size(mesh.x)));
W0 = zeros(size(mesh.x));
W1 = zeros(size(mesh.x));
W2 = zeros(size(mesh.x));
W25 = zeros(size(mesh.x));
W50 = zeros(size(mesh.x));
for i = 1:size(mesh.x,1)
    for j = 1:size(mesh.x,2)
        %Methods from Stelling and Zljlema
         W0(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),0}))*(obj.Eta0(i,j));
         W1(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),1}))*(obj.Eta1(i,j));
         W2(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),2}))*(obj.Eta2(i,j));       
        
%         W0(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),0}))*(obj.Eta0(i,j)+h);
%         W10(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),10}))*(obj.Eta0(i,j)+h);
%         W20(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),20}))*(obj.Eta20(i,j)+h);
%         W25(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),25}))*(obj.Eta25(i,j)+h);
%         W50(i,j) = double(subs(difU,{x,t},{mesh.x(i,j),50}))*(obj.Eta50(i,j)+h);
    end
end
obj.W0 = -W0;
obj.W1 = -W1;
obj.W2 = -W2;
% obj.W25 = -W25;
% obj.W50 = -W50;
%% ���ű��������������
% % obj.W0 = -double(subs(difU,{x,t},{mesh.x,0})).*(obj.Eta0+h);
% % obj.W25 = -double(subs(difU,{x,t},{mesh.x,25})).*(obj.Eta25+h);
% % obj.W50 = -double(subs(difU,{x,t},{mesh.x,50})).*(obj.Eta50+h);
end