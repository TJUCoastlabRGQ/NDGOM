% format long;
% a = 2.5 - sqrt(9.8);
% b = 3.86398-sqrt(9.8*0.611753);
% 
% % h = 1 - 0.4*abs(a)/(abs(a)+abs(b));
% 
% u = 2.5 + (3.86398 - 2.5 )*abs(a)/(abs(a)+abs(b));
% 
% u1 = 1/3*(2.5 + 2*sqrt(9.8));
% 
% format short;

rx = [3 5 7]';
h = [1 2.4 5.1]';
Dr = [0.2 0.4 0.5; 0.37 0.45 0.63; 0.89 0.91 0.75];
q = [3 6 9]';
deltat = 0.93;
rhoreverse = 0.002;
deltat*rhoreverse*rx.*(Dr*((1./h).*(rx.*(Dr*(0.5*h.*q)))));
deltat*rhoreverse*diag(rx)*(Dr*(diag((1./h))*(diag(rx)*(Dr*(0.5*diag(h)*q)))));
deltat*rhoreverse*diag(rx)*Dr*diag(1./h)*diag(rx)*Dr*0.5*diag(h)*q;
(deltat*rhoreverse*diag(rx)*Dr*diag(1./h)*diag(rx)*Dr*0.5*diag(h))*q;



% need to be verified
0.5*deltat*rhoreverse*rx.^2.*((1./h).*(Dr*(h.*q)));
0.5*deltat*rhoreverse*rx.^2.*(1./h).*(Dr*(h.*q))
% wrong
% deltat*rhoreverse*rx.*((1./h).*(rx.*(Dr*(0.5*h)))).*q;
% deltat*rhoreverse*rx.^2.*(Dr*((1./h).*(Dr*(0.5*h.*q))));
% (deltat*rhoreverse*diag(rx)*diag(rx)*Dr*diag(1./h)*Dr*0.5*diag(h))*q;
1/2*h;
1./h;



Lambda = 20;
x0 = 10;
a = 0.01;
h = 15;
c = sqrt(9.8*Lambda/2/pi*tanh(2*pi*h/Lambda));
T = Lambda/c;
tempx = zeros(10/0.01+1,1);
x = 0: 0.01 : 10;
tempdata = zeros(size(x));
k = 2*pi/Lambda;
figure;
for t = 0: 0.01 : 10
 tempdata = 2*a*cos(k.* x)* cos(2*pi/T*t);
 plot(x,tempdata);
 hold on;
 pause(0.05);
 cla;
end

k = 1;
for i = 0 : 0.5 :600
    x(k) = i;
    eta(k) = 2*(sech(sqrt(3*2/4/1000) * (i - 80) ))^2;
    k = k+1;
end
plot(x , eta);

