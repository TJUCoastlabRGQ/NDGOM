function data
close all;
ydata = [1 3 5 7 9];
xdata = [-2 0 2];
for i = 1:5

Title = strcat('WaveTransformOverAnEllipticalShoal\','section',num2str(i),'.fig');
h = open(Title); % open figure
[ measuredX, simulatedX, measuredY, simulatedY ]=getdata(gca);
close(h);
newh = plotdata(measuredX, simulatedX, measuredY, simulatedY);
str = strcat('$$','\rm{Section }',num2str(i), ' (','y=',num2str(ydata(i)),'\rm{m}',')','$$');

% legend({str,'$i$','(','$y=\;ydata(i)\rm{(m)}$'}, 'FontSize',15,'Interpreter','latex');
% legend('$str,i,(,y=\;ydata(i)\rm{(m)}$', 'FontSize',15,'Interpreter','latex');
% legend(str,'Interpreter','latex');
t = text ( 3, 2.2, str,'Interpreter','latex','FontSize',15);
t.FontSize = 15;
end

for i = 1:3

Title = strcat('WaveTransformOverAnEllipticalShoal\','section',num2str(i+5),'.fig');
h = open(Title); % open figure
[ measuredX, simulatedX, measuredY, simulatedY ]=getdata(gca);
close(h);
newh = plotdata(measuredX, simulatedX, measuredY, simulatedY);
str = strcat('$$','\rm{Section }',num2str(i+5), ' (','x=',num2str(xdata(i)),'\rm{m}',')','$$');

% legend({str,'$i$','(','$y=\;ydata(i)\rm{(m)}$'}, 'FontSize',15,'Interpreter','latex');
% legend('$str,i,(,y=\;ydata(i)\rm{(m)}$', 'FontSize',15,'Interpreter','latex');
% legend(str,'Interpreter','latex');
if i~=1
t = text ( 9.6, 2.2, str,'Interpreter','latex','FontSize',15);
else
t = text ( 9.4, 2.2, str,'Interpreter','latex','FontSize',15);
end
t.FontSize = 15;
end
end

function [ measuredX, simulatedX, measuredY, simulatedY ]=getdata(Currentaxis)
D1=get(Currentaxis,'Children'); %get the handle of the line object
XData=get(D1,'XData'); %get the x data
YData=get(D1,'YData'); %get the y data
sizedata = zeros(1,2);
sizedata(1) = numel(XData{1});
sizedata(2) = numel(XData{2});
if sizedata(1)>sizedata(2)
    measuredX = XData{2};
    measuredY = YData{2};
    simulatedX = XData{1};
    simulatedY = YData{1};
else
    measuredX = XData{1};
    measuredY = YData{1};
    simulatedX = XData{2};
    simulatedY = YData{2};    
end
end

function h = plotdata(measuredX, simulatedX, measuredY, simulatedY)
h = figure;
set(gcf,'position',[50,50,1050,400]);
plot(simulatedX, simulatedY,'k','LineWidth',1.5);
hold on;
plot(measuredX, measuredY,'ro','markersize',5);
set(gca,'YLIM',[0, 2.5],'Fontsize',15);
set(gca,'Fontsize',15);
xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
ylabel({'$h/h0$'},'Interpreter','latex');

end