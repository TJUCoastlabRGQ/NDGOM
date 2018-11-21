function CheckProfileResult(obj)
y = [1, 3, 5, 7, 9];
Nout = obj.outputFile.outputStep;
fontsize = 15;
% markersize = 6;

for i = 1:numel(y)
    str = strcat('Section',num2str(i));
    csv = strcat('profile',num2str(i),'.csv');
    filename = [ fileparts( mfilename('fullpath') ), ...
    '/data/',csv];
    data = xlsread(filename);
    figure;
    xpf = -5:0.05:5;
    ypf = y(i).*ones(size(xpf));
    analys = Analysis2d( obj, xpf, ypf );
    temp = analys.InterpOutputResult2GaugePoint( ceil(Nout/1) );
    gaugeResult = temp( :, 1 );
    Initialdata = analys.InterpOutputResult2GaugePoint( 1 );
%     plot(xpf,(gaugeResult./Initialdata(:,1))','LineWidth',3);
    plot(xpf,gaugeResult./obj.amplitude,'LineWidth',3);
    hold on;
    plot(data(:,1),data(:,2),'LineWidth',3);
    legend(str);
end
end