function Solitarywave(obj, mesh)
syms x z t;
rho = 1000;
h = obj.InitialDepth;
a = 0.2 * h;
epsilon = a/h;
alpha = sqrt(3/4*epsilon)*(1-5/8*epsilon +71/128*epsilon^2);
s = sech(alpha*(x*h-100));
t = tanh(alpha*(x*h-100));
eta = h*(1 + epsilon * s^2 - 3/4*epsilon^2*s^2*t^2 + epsilon^3*(5/8*s^2*t^2-101/80*s^4*t^2));
y = eta * (1+z);
u = (1 + 1/2*epsilon - 3/20*epsilon^2 + 3/56*epsilon^3 - epsilon*s^2 + epsilon^2*(-1/4*s^2+s^4+y^2*(3/2*s^2-9/4*s^4)) + ...
    epsilon^3*(19/40*s^2 + 1/5*s^4-6/5*s^6 + y^2*(-3/2*s^2-15/4*s^4+15/2*s^6) + y^4*(-3/8*s^2+45/16*s^4-45/16*s^6)) )*sqrt(obj.gra*h*(1));
w = sqrt(3*epsilon)*y*t*(-1*epsilon*s^2+epsilon^2*(3/8*s^2+2*s^4+y^2*(1/2*s^2-3/2*s^4)) + ...
    epsilon^3*(49/640*s^2 - 17/20*s^4 - 18/5*s^6 + y^2*(-13/16*s^2-25/16*s^4+15/2*s^6)+y^4*(-3/40*s^2+9/8*s^4-27/16*s^6)))*sqrt(obj.gra*h*(1));
p = (1-y + epsilon*s^2 + epsilon^2*(3/4*s^2-3/2*s^4+y^2*(-3/2*s^2+9/4*s^4)) + epsilon^3*(-1/2*s^2-19/20*s^4+11/5*s^6+y^2*(3/4*s^2+39/8*s^4-33/4*s^6)+y^4*(3/8*s^2-45/16*s^4+45/16*s^6)))*rho*obj.gra*h;

t = 0;
x = mesh.mesh2d.x;
obj.H0 = eval(eta);
z = -0.5;
obj.CEU0 = eval(u);
obj.CEP0 = eval(p);
obj.CEW0 = eval(w);

x = mesh.x;
z = mesh.z;
obj.U0 = eval(u);
obj.W0 = eval(w);
% 
% 
% % t = 10;
% % x = mesh.mesh2d.x + 10*sqrt(obj.gra*h*(1+epsilon));
% x = mesh.mesh2d.x + 10*sqrt(obj.gra*h*(1));
% obj.H10 = eval(eta);
% z = -0.5;
% obj.CEU10 = eval(u);
% obj.CEP10 = eval(p);
% obj.CEW10 = eval(w);
% 
% % t = 20;
% x = mesh.mesh2d.x + 20*sqrt(obj.gra*h*(1));
% obj.H20 = eval(eta);
% z = -0.5;
% obj.CEU20 = eval(u);
% obj.CEP20 = eval(p);
% obj.CEW20 = eval(w);
% 
% 
% % t = 30;
% x = mesh.mesh2d.x + 30*sqrt(obj.gra*h*(1));
% obj.H30 = eval(eta);
% z = -0.5;
% obj.CEU30 = eval(u);
% obj.CEP30 = eval(p);
% obj.CEW30 = eval(w);
% 
% 
% % syms x z t;
% % h = obj.InitialDepth;
% % a = h*0.2;
% % x0 = 1/2*obj.ChLength;
% % phi = 1/2/h*sqrt(3*a/(h+a));
% % c = sqrt(obj.gra*(h+a));
% % 
% % eta = 4*a*exp(-2*phi*(x-x0-c*t))/(1+exp(-2*phi*(x-x0-c*t)))^2;
% % % v = c*(eta-h)/eta;
% % v = eta * c/(eta + h);
% % w = -1*diff(v,x)*(eta + h)*(1 + z);
% % 
% % xmax = max( max(mesh.mesh2d.x) );
% % xmin = min( min(mesh.mesh2d.x) );
% % T = mod(10, (xmax - xmin)/c);
% % Index = mesh.mesh2d.x - xmin >= c * T;
% % obj.H10(Index) = double(subs(eta,{x,t},{ mesh.mesh2d.x(Index) - c * T,  0})) + h;
% % Index = mesh.mesh2d.x - xmin <= c * T;
% % obj.H10(Index) = double(subs(eta,{x,t},{ mesh.mesh2d.x(Index) - c * T + xmax - xmin,  0})) + h;
% % 
% % Index = mesh.x - xmin >= c * T;
% % obj.U10(Index) = double(subs(v,{x,t},{ mesh.x(Index) - c * T,  0}));
% % % obj.P5(Index) = double(subs(P,{x,t},{ mesh.x(Index) - c * T,  0}));
% % obj.W10(Index) = double(subs(w,{x,t,z},{ mesh.x(Index) - c * T,  0, mesh.z(Index)}));
% % Index = mesh.x - xmin <= c * T;
% % obj.U10(Index) = double(subs(v,{x,t},{ mesh.x(Index) - c * T + xmax - xmin,  0}));
% % % obj.P5(Index) = double(subs(P,{x,t},{ mesh.x(Index) - c * T,  0}));
% % obj.W10(Index) = double(subs(w,{x,t,z},{ mesh.x(Index) - c * T + xmax - xmin,  0, mesh.z(Index)}));
% % 
% % 
% % T = mod(20, (xmax - xmin)/c);
% % Index = mesh.mesh2d.x - xmin >= c * T;
% % obj.H20(Index) = double(subs(eta,{x,t},{ mesh.mesh2d.x(Index) - c * T,  0})) + h;
% % Index = mesh.mesh2d.x - xmin <= c * T;
% % obj.H20(Index) = double(subs(eta,{x,t},{ mesh.mesh2d.x(Index) - c * T + xmax - xmin,  0})) + h;
% % 
% % Index = mesh.x - xmin >= c * T;
% % obj.U20(Index) = double(subs(v,{x,t},{ mesh.x(Index) - c * T,  0}));
% % % obj.P5(Index) = double(subs(P,{x,t},{ mesh.x(Index) - c * T,  0}));
% % obj.W20(Index) = double(subs(w,{x,t,z},{ mesh.x(Index) - c * T,  0, mesh.z(Index)}));
% % Index = mesh.x - xmin <= c * T;
% % obj.U20(Index) = double(subs(v,{x,t},{ mesh.x(Index) - c * T + xmax - xmin,  0}));
% % % obj.P5(Index) = double(subs(P,{x,t},{ mesh.x(Index) - c * T,  0}));
% % obj.W20(Index) = double(subs(w,{x,t,z},{ mesh.x(Index) - c * T + xmax - xmin, 0, mesh.z(Index)}));
% % 
% % T = mod(30, (xmax - xmin)/c);
% % Index = mesh.mesh2d.x - xmin >= c * T;
% % obj.H30(Index) = double(subs(eta,{x,t},{ mesh.mesh2d.x(Index) - c * T,  0})) + h;
% % Index = mesh.mesh2d.x - xmin <= c * T;
% % obj.H30(Index) = double(subs(eta,{x,t},{ mesh.mesh2d.x(Index) - c * T + xmax - xmin,  0})) + h;
% % 
% % Index = mesh.x - xmin >= c * T;
% % obj.U30(Index) = double(subs(v,{x,t},{ mesh.x(Index) - c * T,  0}));
% % % obj.P5(Index) = double(subs(P,{x,t},{ mesh.x(Index) - c * T,  0}));
% % obj.W30(Index) = double(subs(w,{x,t,z},{ mesh.x(Index) - c * T,  0, mesh.z(Index)}));
% % Index = mesh.x - xmin <= c * T;
% % obj.U30(Index) = double(subs(v,{x,t},{ mesh.x(Index) - c * T + xmax - xmin,  0}));
% % % obj.P5(Index) = double(subs(P,{x,t},{ mesh.x(Index) - c * T,  0}));
% % obj.W30(Index) = double(subs(w,{x,t,z},{ mesh.x(Index) - c * T + xmax - xmin,  0, mesh.z(Index)}));
% % 
% % x = mesh.mesh2d.x;
% % t = 0;
% % obj.H0 = eval(eta) + h;
% % 
% % x = mesh.x;
% % z = mesh.z;
% % obj.U0 = eval(v);
% % obj.W0 = eval(w);

% syms x z t;
% rho = 1000;
% h = obj.InitialDepth;
% a = 0.2 * h;
% epsilon = a/h;
% alpha = sqrt(3/4*epsilon)*(1-5/8*epsilon +71/128*epsilon^2);
% C = sqrt(obj.gra * h)*(1 + epsilon - 1/20 * epsilon^2 - 3/70*epsilon^3)^(1/2);
% s = sech(alpha*(x-30 - C*t)/h);
% q = tanh(alpha*(x-30 - C*t)/h);
% 
% eta = h*(epsilon * s^2 - 3/4*epsilon^2*s^2*q^2 + epsilon^3*(5/8*s^2*q^2-101/80*s^4*q^2));
% y = (eta + h) * (1+z);
% 
% u = sqrt(obj.gra * h) * ( epsilon * s^2 + epsilon^2 * (-3/4*s^2 + s^2*q^2 + (y/h)^2*(3/4*s^2-9/4*s^2*q^2)) + ...
%     epsilon^3 * (21/40*s^2 - s^2*q^2 - 6/5*s^4*q^2 + (y/h)^2*(-9/4*s^2+15/4*s^2*q^2+15/2*s^4*q^2) + (y/h)^4*(3/8*s^2 - 45/16*s^4*q^2)) );
% 
% w = sqrt(3*epsilon)*(y/h)*s^2*q*(epsilon + epsilon^2*(-3/8 - 2*s^2 + (y/h)^2*(-1/2 + 3/2*s^2)) + epsilon^3 * ...
%     (-49/640 - 17/20*s^2 - 18/5*s^4 + (y/h)^2*(-13/16 - 25/16 * s^2 + 15/2 * s^4) + (y/h)^4 * (-3/40 + 9/8*s^2 - 27/16*s^4))) * sqrt(obj.gra * h);
% 
% 
% t = 0;
% x = mesh.mesh2d.x;
% obj.H0 = eval(eta) + h;
% z = -0.5;
% obj.CEU0 = eval(u);
% % obj.CEP0 = eval(p);
% obj.CEW0 = eval(w);
% 
% x = mesh.x;
% z = mesh.z;
% obj.U0 = eval(u);
% obj.W0 = eval(w);
% % Index = (x>=150+7.715 | x<=150-7.715);
% % obj.U0(Index) = 0;
% % obj.W0(Index) = 0;
% 
% % t = 0.5;
% % z = -0.5;
% % x = mesh.mesh2d.x;
% % obj.CEW05 = eval(w);
% 
% t = 10;
% % x = mesh.mesh2d.x + 10*sqrt(obj.gra*h*(1+epsilon));
% x = mesh.mesh2d.x;
% obj.H10 = eval(eta) + h;
% z = -0.5;
% obj.CEU10 = eval(u);
% % obj.CEP10 = eval(p);
% obj.CEW10 = eval(w);
% 
% t = 20;
% x = mesh.mesh2d.x;
% obj.H20 = eval(eta) + h;
% z = -0.5;
% obj.CEU20 = eval(u);
% % obj.CEP20 = eval(p);
% obj.CEW20 = eval(w);
% 
% t = 30;
% x = mesh.mesh2d.x;
% obj.H30 = eval(eta) + h;
% z = -0.5;
% obj.CEU30 = eval(u);
% % obj.CEP30 = eval(p);
% obj.CEW30 = eval(w);
end