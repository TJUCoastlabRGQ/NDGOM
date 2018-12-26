function TsunamiRunup2dPostprocess
clear,clc
close all
load NhydroSolver;
time = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 22.5];

PostProcess = NdgPostProcess(Solver.meshUnion(1),strcat('NhydroTsunamiRunup2d','.',num2str(Solver.Nmesh),'-','1','/','TsunamiRunup2d'));
%% draw 3D view
% draw3Dview(PostProcess, time, Solver);

%% plot the topography
drawTopography(Solver);

%% plot the boundary condition
drawBoundaryCondition(Solver);

PointTimeSeriesProcess

end
%% function to draw the 3D view
function draw3Dview(PostProcess, time, Solver)
outputTime = ncread( PostProcess.outputFile{1}, 'time' );
for i  = 1:numel(time)
    figure;
    [~, index] = sort(abs(outputTime - time(i)));
    Visual = makeVisualizationFromNdgPhys( Solver );
    fphys = PostProcess.accessOutputResultAtStepNum(index(1));
%     set(gca,'XLim',[3,5]);
%     set(gca,'YLim',[0.5,3.5]);
    Visual.drawResult(fphys{1}(:, :, 1)+Solver.fphys{1}(:, :, 4) - Solver.hmin);
    set(Visual.drawHandle,'edgecolor','none','facecolor','interp');
    light('Position',[3 1.5 1.5],'Style','local');
    
    Visual.drawResult(Solver.fphys{1}(:, :, 4));
    set(Visual.drawHandle,'facecolor',[0.5,0.5,0.5],'edgecolor','none');
    set(gcf, 'position',[100 100 1200 830]);
    
    view(240,40);
    set(gca,'Fontsize',18);
    xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex','Fontsize',20);
    ylabel({'$y\;\rm{(m)}$'},'Interpreter','latex','Fontsize',20);
    zlabel({'$h+b\;\rm{(m)}$'},'Interpreter','latex','Fontsize',20);
    h = colorbar('southoutside');
    h.Label.FontSize = 20;
    h.Limits = [-0.1 0.1];
    h.Position = [0.075 0.035 0.85 0.02];
end
end
%% function to draw the topography data
function drawTopography(Solver)

fp = fopen([ 'D:\PhdResearch\Application\SWE\SWE2d\Benchmark\@TsunamiRunup2d', ...
    '/mesh/Benchmark_2_Bathymetry.txt']);
fgetl(fp);
data = fscanf(fp, '%e %e %e', [3, inf]);
fclose(fp);

x = linspace(min(min(Solver.meshUnion(1).x)),max(max(Solver.meshUnion(1).x)));
y = linspace(min(min(Solver.meshUnion(1).y)),max(max(Solver.meshUnion(1).y)));
[X,Y] = meshgrid(x,y);

interp = scatteredInterpolant( ...
    data(1,:)', data(2,:)', -data(3,:)','linear');
bot = interp( X, Y );
contour(X,Y,bot,25,'Linewidth',1.5);
handle = colorbar('southoutside');
handle.Label.FontSize = 20;
h.Position = [0.005 0.005 0.85 0.01];
set(gca,'Fontsize',15);
xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex','Fontsize',20);
ylabel({'$y\;\rm{(m)}$'},'Interpreter','latex','Fontsize',20);
set(gcf, 'position',[100 100 1000 600]);
hold on;

PositionX = [4.521 4.521 4.521];
PositionY = [1.196 1.696 2.196];
Str = {'$G5$','$G7$','$G9$'};
plot(PositionX, PositionY, 'o', 'markersize',10, 'markerfacecolor',[0 0 0], 'markeredgecolor','k');
for i = 1:numel(PositionX)
    text(PositionX(i) - 0.3, PositionY(i), Str{i}, 'Interpreter','latex','fontsize', 15);
end
end
%% function to plot boundary conditon
function drawBoundaryCondition(Solver)
figure;
plot(Solver.boundaryInfo(1).time,Solver.boundaryInfo(1).extElevation,'Linewidth',1.5);
set(gcf, 'position',[100 100 1000 600]);
set(gca,'Fontsize',15);
xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex','Fontsize',20);
ylabel({'$h+b\;\rm{(m)}$'},'Interpreter','latex','Fontsize',20);
grid on;
end
% %% function to postprocess the point time series data

function PointTimeSeriesProcess
%% function to postprocess the point time series data
close all
str = 'hydrostatic-rk45.fig';
str1 = 'Non-hydrostatic-rk45.fig';
h = openfig(str,'reuse'); % open figure
h1 = openfig(str1,'reuse'); % open figure 

num = numel(h.Children);
Xposition = zeros(1,num);
Yposition = zeros(1,num);
%% read coordinate data
for i = 1:num
Xposition(i) = h.Children(i).Position(1);    
Yposition(i) = h.Children(i).Position(2);    
end
%% find the x postion of the data figure
flag = 1;
i = 1;
while flag
tempindex = find(abs(Xposition - Xposition(i))<0.001);
Ydata = zeros(size(tempindex));
for k = 1:numel(tempindex)
    Ydata(k) = Yposition(tempindex(k));
end
 if numel(tempindex)>1
     flag = 0;
 end
 i = i+1;
end

[index, ~] = sort(tempindex,'descend');

%% get data
[HXdata, HYdata,MeasuredXdata,MeasuredYdata] = getpointdata(h, index);
[NhXdata, NhYdata,~,~] = getpointdata(h1, index);

%% plot data
for i = 1:numel(index)
    figure;
    set(gcf, 'position',[100 100 1000 400]);
    plot(MeasuredXdata{i},MeasuredYdata{i},'r-',...
    HXdata{i},HYdata{i},'k',...
    NhXdata{i},NhYdata{i},'g','Linewidth',1.5);
    set(gca,'Xlim',[0 22.5]);
    legend({'Measured','Hydro','Nonhydro'},'position',[0 0.75 0.4 0.15]);
    legend('boxoff')
    xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex','Fontsize',15);
    ylabel({'$h+b\;\rm{(m)}$'},'Interpreter','latex','Fontsize',15);
    set(gca,'Fontsize',15);
    
    h1=axes('position',[0.4 0.6 0.2 0.25]);
    axis(h1);
    plot(MeasuredXdata{i},MeasuredYdata{i},'r-',...
    HXdata{i},HYdata{i},'k',...
    NhXdata{i},NhYdata{i},'g','Linewidth',1.5);

%     set(gca,'Xlim',[17 20],'Ylim',[0.02 0.04],'Fontsize',15);  %FOR G5
%     set(gca,'Xlim',[16.5 18.5],'Ylim',[0.027 0.04],'Fontsize',15);  %FOR G7
    set(gca,'Xlim',[16.5 18.5],'Ylim',[0.032 0.045],'Fontsize',15);  %FOR G9
end
end
%%
function [Xdata, Ydata, MeasuredX, MeasuredY] = getpointdata(handle, index)
Xdata = cell(numel(index),1);MeasuredX = cell(numel(index),1);
Ydata = cell(numel(index),1);MeasuredY = cell(numel(index),1);
for i = 1:numel(index)
ca = get(handle.Children(index(i)),'Children');
tempXData =get(ca,'XData');
tempYData =get(ca,'YData');
datasize = zeros(size(tempXData));
  for j = 1:numel(tempXData)
      datasize(j) = numel(tempXData{j});
  end
  [~, ascendindex] = sort(datasize,'ascend');
  Xdata{i} = tempXData{ascendindex(1)};
  MeasuredX{i} = tempXData{ascendindex(2)};
  Ydata{i} = tempYData{ascendindex(1)};
  MeasuredY{i} = tempYData{ascendindex(2)};  
end
%   for j = 1:numel(tempXData)
%      Xdata{j} = tempXData{ascendindex(j)};
%      Xdata{j+numel(index)} = tempXData{ascendindex(j)};
%      Ydata{j} = tempYData{ascendindex(j)};
%   end
end
