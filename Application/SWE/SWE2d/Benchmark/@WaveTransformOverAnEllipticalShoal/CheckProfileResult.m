function CheckProfileResult(obj)
y = [1, 3, 5, 7, 9];
% x = [-2 0 2];
Nout = obj.outputFile.outputStep;
% fontsize = 15;
% markersize = 6;

for i = 1:numel(y)
% for i = 1:numel(x)
%     str = strcat('Section',num2str(i+5));
%     csv = strcat('profile',num2str(i+5),'.csv');
    str = strcat('Section',num2str(i));
    csv = strcat('profile',num2str(i),'.csv');
    filename = [ fileparts( mfilename('fullpath') ), ...
    '/data/',csv];
    data = xlsread(filename);
    figure;
    xpf = -5:0.05:5;
    ypf = y(i).*ones(size(xpf));
%     ypf = 0:0.05:12;
%     xpf = x(i).*ones(size(ypf));
    analys = Analysis2d( obj, xpf,  ypf);
    
    tempgaugeResult = zeros(numel(xpf),ceil(Nout/2)-ceil(Nout/4)+1);
    k = 1;
    for t = ceil(Nout/4):ceil(1*Nout/2)
       tempdata = analys.InterpOutputResult2GaugePoint( t ); 
       tempgaugeResult( :,k ) = tempdata( :,1 );
       k = k+1;
    end  
    depthdata = sort( tempgaugeResult, 2 );
    gaugeResult = depthdata(:,end) - depthdata(:,1);
    
%     Initialdata = analys.InterpOutputResult2GaugePoint( 1 );
%     plot(xpf,(gaugeResult./Initialdata(:,1))','LineWidth',3);
    plot(xpf,gaugeResult./obj.amplitude/2,'LineWidth',3);
%     plot(ypf,gaugeResult./obj.amplitude/2,'LineWidth',3);
    hold on;
    plot(data(:,1),data(:,2),'LineWidth',3);
    legend(str);
end
end