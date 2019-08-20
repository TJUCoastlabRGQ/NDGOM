function drawGaugeResult( obj )

% gx = [ 6.71 9.36 10.36 12.96 15.56]; %case 0.096
% gx = [ 7.45 9.36 10.36 12.96 15.56]; %case 0.181
gx = [ 5.91 9.36 10.36 12.96 15.56]; %case 0.045
gy = [ 13.05 13.8 13.8 11.22 13.8];
GaugePoint = [3 6 9 16 22];

conicalPos = NdgPostProcess(obj.meshUnion(1),strcat('SolitaryWaveRunUpConicalIsland','/','SolitaryWaveRunUpConicalIsland'));
[ result ] = conicalPos.interpolateOutputResultToGaugePoint( gx, gy, gx );
[ time ] = ncread( conicalPos.outputFile{1}, 'time' );
[ gaugeValue ] = conicalPos.interpolatePhysFieldToGaugePoint( ...
    obj.fphys, gx, gy, gx );
[ bot ] = gaugeValue(:, 4)';
[ path, name, ext ] = fileparts( mfilename('fullpath') );

% for paper

for i = 1:numel(gx)
    str = strcat('\result\case 0.045\Gauge',num2str(GaugePoint(i)),'.csv');
    pathstr = strcat(path,str);
    titlestr = strcat('Gauge ',32 ,num2str(GaugePoint(i)));
    data = xlsread(pathstr);
    figure;
    set(gcf,'position',[50,50,265,150]);
%     plot(time + 22,( result(:, 1, i) - 0.32 + bot(i)     )*100,'k','LineWidth',1.5); %case 0.096
    plot(time + 25.6,( result(:, 1, i) - 0.32 + bot(i) )*100,'k','LineWidth',1.5); %case 0.181
    hold on;
    plot(data(:,1),( data(:,2) )*100,'ro','markersize',3);
    ylim([-5,5]);
    ylabel({'$\eta \;\rm {(cm)}$'},'Interpreter','latex');
    xlabel({'$t \;\rm {(s)}$'},'Interpreter','latex');
    set(gca,'Fontsize',10);
    title(titlestr,'position',[40,2.8],'Fontsize',10);
    xlim([ 25, 45 ]);
end
% for thesis

for i = 1:numel(gx)
    str = strcat('\result\case 0.045\Gauge',num2str(GaugePoint(i)),'.csv');
    pathstr = strcat(path,str);
    titlestr = strcat('Gauge',32,num2str(GaugePoint(i)));
    data = xlsread(pathstr);
    figure;
    set(gcf,'position',[50,50,1000,400]);
%     plot(time + 22,( result(:, 1, i) - 0.32 + bot(i) )*100,'k','LineWidth',2.5);%case 0.096
%     plot(time + 21.2,( result(:, 1, i) - 0.32 + bot(i) )*100,'k','LineWidth',2.5);%case 0.181
     plot(time + 25.6,( result(:, 1, i) - 0.32 + bot(i) )*100,'k','LineWidth',2.5);%case 0.045
    hold on;
    plot(data(:,1),( data(:,2) )*100,'ro','markersize',6.5);
    ylim([-5,5]);
    ylabel({'$\eta \;\rm {(cm)}$'},'Interpreter','latex');
    xlabel({'$t \;\rm {(s)}$'},'Interpreter','latex');
    set(gca,'Fontsize',18);
    title(titlestr,'position',[42,3.5],'Fontsize',18);%case 0.045
    xlim([ 25, 45 ]);
end


end
