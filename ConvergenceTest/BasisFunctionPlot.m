function BasisFunctionPlot
%BASISFUNCTIONPLOT 此处显示有关此函数的摘要
%   此处显示详细说明
close all;
clear,clc;
Solver = EkmanSpiral(2,2,2000000,2);
for i = 1:size(Solver.mesh2d.cell.V,1)
    figure;
    surf(Solver.mesh2d.cell.Drawr, Solver.mesh2d.cell.Draws, reshape(Solver.mesh2d.cell.DrawV(:,i),size(Solver.mesh2d.cell.Drawr,1),size(Solver.mesh2d.cell.Drawr,2)),...
        'Facecolor','white');
    str = ['$\varphi_',num2str(i),'(r,s)$'];
    title(str, 'interpreter','latex','Fontsize',20);
    set(gca,'fontsize',20,'linewidth',1.2);
    zlim([-2,4]);
    grid on
end
close all;
r = [-1 0 1 -1 0 -1]';
s = [-1 -1 -1 0 0 1]';
% For Lagrange interpolation function
Data = (Solver.mesh2d.cell.V)'\Solver.mesh2d.cell.DrawV';
for i = 1:size(Solver.mesh2d.cell.V,1)
    figure;
    surf(Solver.mesh2d.cell.Drawr, Solver.mesh2d.cell.Draws, reshape(Data(i,:),size(Solver.mesh2d.cell.Drawr,1),size(Solver.mesh2d.cell.Drawr,2)),...
        'Facecolor','white');
    hold on;
    z = zeros(6,1);
    z(i) = 1;
    scatter3(r,s,z,80,'black','filled');
    str = ['$l_',num2str(i),'(r,s)$'];
    title(str, 'interpreter','latex','Fontsize',20);
    set(gca,'fontsize',20,'linewidth',1.2);
    zlim([-0.5,1.5]);
    grid on
end
close all
%For mode function
for i = 1:size(Solver.mesh2d.InnerEdge.cell.V,1)
    figure;
    plot(Solver.mesh2d.InnerEdge.cell.Drawr, Solver.mesh2d.InnerEdge.cell.DrawV(:,i)','black','Linewidth',1.5);
    ylim([-1.5,2]);
    str = ['$\varphi_',num2str(i),'(r)$'];
    title(str, 'interpreter','latex','Fontsize',20);
    set(gca,'fontsize',20,'linewidth',1.2);
end
close all
% For Lagrange interpolation function
r = [-1 0 1]';
Data = (Solver.mesh2d.InnerEdge.cell.V)'\Solver.mesh2d.InnerEdge.cell.DrawV';
for i = 1:size(Solver.mesh2d.InnerEdge.cell.V,1)
    figure;
    plot(Solver.mesh2d.InnerEdge.cell.Drawr, Data(i,:),'black','Linewidth',1.5);

    hold on;
    z = zeros(3,1);
    z(i) = 1;
    scatter(r,z,80,'black','filled');
    str = ['$l_',num2str(i),'(r)$'];
    title(str, 'interpreter','latex','Fontsize',20);
    set(gca,'fontsize',20,'linewidth',1.2);
    ylim([-0.5,1.5]);
    grid on
end
close all
end

