function CheckProfileResult(obj)
y = [1, 3, 5, 7, 9]; x = [-2 0 2];
% fontsize = 15;
% markersize = 6;

YprofileProcess(obj, y);
XprofileProcess(obj, x);

end

function  YprofileProcess(obj, y)
% Nout = obj.outputFile.outputStep;
WaveTransformPos = NdgPostProcess(obj.meshUnion(1).mesh2d,strcat('WaveTransformOverAnEllipticalShoal','/','WaveTransformOverAnEllipticalShoal'));
[ time ] = ncread( WaveTransformPos.outputFile{1}, 'time' );
Nout = numel(time);

for i = 1:numel(y)
    str = strcat('$Section',32,num2str(i),'$');
    csv = strcat('profile',num2str(i),'.csv');
    filename = [ fileparts( mfilename('fullpath') ), ...
    '/data/',csv];
    data = xlsread(filename);
  
    xpf = -5:0.1:5;
    ypf = y(i).*ones(size(xpf));
    
    [ result ] = WaveTransformPos.interpolateOutputResultToGaugePoint( xpf, ypf, xpf );
    index = ceil(25* Nout/27.94)+1:ceil(27.94* Nout/27.94);
    
%     analys = Analysis2d( obj, xpf,  ypf);
    
%     tempgaugeResult = zeros(numel(xpf),ceil(29* Nout/35)-ceil(33* Nout/35)+1);
    tempgaugeResult = zeros(numel(xpf),numel(index));
%     k = 1;
%     for t = ceil(30* Nout/35):ceil(35* Nout/35)
%        tempdata = analys.InterpOutputResult2GaugePoint( t ); 
%        tempgaugeResult( :,k ) = tempdata( :,1 );
%        k = k+1;
%     end  
    for j = 1:numel(xpf)
        tempgaugeResult(j,:) = result(index, 1, j);
    end
    
    depthdata = sort( tempgaugeResult, 2 );
    gaugeResult = depthdata(:,end) - depthdata(:,1);
    
    
    figure;
    set(gcf,'position',[50,50,800,300]);
    plot(xpf,gaugeResult./obj.amplitude/2,'LineWidth',1.5);
    hold on;
    plot(data(:,1),data(:,2),'ro','markersize',5);
    xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
    ylabel({'$h/H_0$'},'Interpreter','latex');
    box on;
    legend({'$Simulated$','$Measured$'},'Interpreter','Latex');
    legend('boxoff')
    text(-4,2,str,'Fontsize',14,'Interpreter','latex');
    xlim([-5,5]);
    ylim([0,2.5]);
    set(gca,'Fontsize',12);
end
end

function  XprofileProcess(obj, x)
% Nout = obj.outputFile.outputStep;
% Nout = 1652;
% Nout = 1616;
WaveTransformPos = NdgPostProcess(obj.meshUnion(1).mesh2d,strcat('WaveTransformOverAnEllipticalShoal','/','WaveTransformOverAnEllipticalShoal'));
[ time ] = ncread( WaveTransformPos.outputFile{1}, 'time' );
Nout = numel(time);

for i = 1:numel(x)
    str = strcat('$Section',32,num2str(i+5),'$');
    csv = strcat('profile',num2str(i+5),'.csv');
    filename = [ fileparts( mfilename('fullpath') ), ...
    '/data/',csv];
    data = xlsread(filename);
  
    ypf = 0:0.1:12;
    xpf = x(i).*ones(size(ypf));
    
    [ result ] = WaveTransformPos.interpolateOutputResultToGaugePoint( xpf, ypf, xpf );
    index = ceil(25* Nout/27.94)+1:ceil(27.94* Nout/27.94);  
    
%     analys = Analysis2d( obj, xpf,  ypf);
    
    tempgaugeResult = zeros(numel(xpf),numel(index));
    
%     k = 1;
%     for t = ceil(30* Nout/35):ceil(35* Nout/35)
%        tempdata = analys.InterpOutputResult2GaugePoint( t ); 
%        tempgaugeResult( :,k ) = tempdata( :,1 );
%        k = k+1;
%     end  

    for j = 1:numel(ypf)
        tempgaugeResult(j,:) = result(index, 1, j);
    end
    
    depthdata = sort( tempgaugeResult, 2 );
    gaugeResult = depthdata(:,end) - depthdata(:,1);
    
    
    figure;
    set(gcf,'position',[50,50,800,300]);
    hold on;
    plot(ypf,gaugeResult./obj.amplitude/2,'LineWidth',1.5);
    plot(data(:,1),data(:,2),'o','markersize',4.5);
    xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
    ylabel({'$h/H_0$'},'Interpreter','latex');
    legend({'$Simulated$','$Measured$'},'Interpreter','Latex');
    legend('boxoff');
    text(1.2,2,str,'Fontsize',14,'Interpreter','latex');
    xlim([0,12]);
    ylim([0,2.5]);
    set(gca,'Fontsize',12);    
end
end