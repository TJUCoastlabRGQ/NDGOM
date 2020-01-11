function CheckProfileResult(obj)
y = [1, 3, 5, 7, 9]; x = [-2 0 2];
% fontsize = 15;
% markersize = 6;

% YprofileProcess(obj, y);
XprofileProcess(obj, x);

end

function  YprofileProcess(obj, y)
% Nout = obj.outputFile.outputStep;
Nout = 1616;
for i = 1:numel(y)
    str = strcat('Section',num2str(i));
    csv = strcat('profile',num2str(i),'.csv');
    filename = [ fileparts( mfilename('fullpath') ), ...
    '/data/',csv];
    data = xlsread(filename);
  
    xpf = -5:0.05:5;
    ypf = y(i).*ones(size(xpf));
    analys = Analysis2d( obj, xpf,  ypf);
    
%     tempgaugeResult = zeros(numel(xpf),ceil(29* Nout/35)-ceil(33* Nout/35)+1);
    tempgaugeResult = zeros(numel(xpf),ceil(35* Nout/35)-ceil(30* Nout/35)+1);
    k = 1;
    for t = ceil(30* Nout/35):ceil(35* Nout/35)
       tempdata = analys.InterpOutputResult2GaugePoint( t ); 
       tempgaugeResult( :,k ) = tempdata( :,1 );
       k = k+1;
    end  
    depthdata = sort( tempgaugeResult, 2 );
    gaugeResult = depthdata(:,end) - depthdata(:,1);
    
    
    figure;
    set(gcf,'position',[50,50,1050,400]);
    set(gca,'YLIM',[0, 2.5],'Fontsize',12);
    set(gca,'XLIM',[-5, 5],'Fontsize',12);
    plot(xpf,gaugeResult./obj.amplitude/2,'LineWidth',1.5);
    hold on;
    plot(data(:,1),data(:,2),'o','markersize',3);
    xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
    ylabel({'$h/h0$'},'Interpreter','latex');
    legend(str);
    legend('boxoff');
end
end

function  XprofileProcess(obj, x)
% Nout = obj.outputFile.outputStep;
% Nout = 1652;
Nout = 1616;
for i = 1:numel(x)
    str = strcat('Section',num2str(i+5));
    csv = strcat('profile',num2str(i+5),'.csv');
    filename = [ fileparts( mfilename('fullpath') ), ...
    '/data/',csv];
    data = xlsread(filename);
  
    ypf = 0:0.15:12;
    xpf = x(i).*ones(size(ypf));
    analys = Analysis2d( obj, xpf,  ypf);
    
    tempgaugeResult = zeros(numel(xpf),ceil(35* Nout/35)-ceil(30* Nout/35)+1);
    k = 1;
    for t = ceil(30* Nout/35):ceil(35* Nout/35)
       tempdata = analys.InterpOutputResult2GaugePoint( t ); 
       tempgaugeResult( :,k ) = tempdata( :,1 );
       k = k+1;
    end  
    depthdata = sort( tempgaugeResult, 2 );
    gaugeResult = depthdata(:,end) - depthdata(:,1);
    
    
    figure;
    set(gcf,'position',[50,50,1050,400]);
    set(gca,'YLIM',[0, 2.5],'Fontsize',12);
    set(gca,'XLIM',[0, 12],'Fontsize',12);
    hold on;
    plot(ypf,gaugeResult./obj.amplitude/2,'LineWidth',1.5);
    plot(data(:,1),data(:,2),'o','markersize',3);
    xlabel({'$y\;\rm{(m)}$'},'Interpreter','latex');
    ylabel({'$h/h0$'},'Interpreter','latex');
    legend(str);
    legend('boxoff');
end
end