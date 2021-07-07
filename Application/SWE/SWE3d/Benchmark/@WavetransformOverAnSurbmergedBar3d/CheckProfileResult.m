function CheckProfileResult(obj)
y = 0.025; x = [10.5 12.5 13.5 14.5 15.7 17.3 19 21];
d = [ 0.175 0.1 0.1 0.15 0.27 0.4 0.398 0.318];
PostProcess = NdgPostProcess(obj.meshUnion(1).mesh2d,'Result/WavetransformOverAnSurbmergedBar3d/2d/WavetransformOverAnSurbmergedBar3d');
Ntime = PostProcess.Nt;
outputTime = ncread( PostProcess.outputFile{1}, 'time' );

Index = find( 33 < outputTime & outputTime < 39 );

Time = outputTime(Index);
Data = zeros(numel(x), numel(Time));
j = 0;

for i = 1:numel(x)
    j = j + 1;
    t = 1;
    for k = 1:numel(Index)
        tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  x(j), y, 0, Index(k) )-d(j);
        
        Data(i,k) = tempdata(1) * 100;
        t = t+1;
    end
    str = ['Gauge',blanks(1),num2str(i+3)];
    csv = strcat(num2str(i + 3),'.csv');
    filename = [ fileparts( mfilename('fullpath') ), ...
        '/data/',csv];
    data = xlsread(filename);
    
    figure;
    set(gcf,'position',[50,50,1050,400]);
    plot(Time',Data(i,:),'LineWidth',1.5);
    hold on;
    plot(data(:,1),data(:,2),'o','markersize',6, 'LineWidth',1.5);
    set(gca,'YLIM',[-2, 4],'Fontsize',12);
    set(gca,'XLIM',[33, 39],'Fontsize',12);
    xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex');
    ylabel({'$\eta\;\rm{(cm)}$'},'Interpreter','latex');
    set(gca,'Fontsize',15);
    title(str,'position',[33.5,3.2],'Fontsize',15);
end

end