function runNonhydroStaticSolver3dTest
% For this test, the covergence rate means nothing, we only use this
% function to calculate the partial derivative.
M = [4 8 16 32 64];
Mz = [2 4 8 16 32];
Order = [1 2];

Nmesh = numel(M);
Ndeg = numel(Order);
len = zeros(Nmesh, Ndeg);
dofs = zeros(Nmesh, Ndeg);

linewidth = 1.5;
markersize = 8;

linestyle = '-';

color = {'k', 'r'};  %black for order one, red for order two
marker = {'o','+','*','.','x','s','d','^','v','>','<'};%circle for the computed field

PSPXErrInf = zeros(Nmesh, Ndeg); PSPYErrInf = zeros(Nmesh, Ndeg);SPSPXErrInf = zeros(Nmesh, Ndeg);
PSPXErr2 = zeros(Nmesh, Ndeg); PSPYErr2 = zeros(Nmesh, Ndeg);SPSPXErr2 = zeros(Nmesh, Ndeg);
PSPXErr1 = zeros(Nmesh, Ndeg); PSPYErr1 = zeros(Nmesh, Ndeg);SPSPXErr1 = zeros(Nmesh, Ndeg);
SPSPYErrInf = zeros(Nmesh, Ndeg);SQPSPXErrInf = zeros(Nmesh, Ndeg);SQPSPYErrInf = zeros(Nmesh, Ndeg);
SPSPYErr2 = zeros(Nmesh, Ndeg);SQPSPXErr2 = zeros(Nmesh, Ndeg);SQPSPYErr2 = zeros(Nmesh, Ndeg);
SPSPYErr1 = zeros(Nmesh, Ndeg);SQPSPXErr1 = zeros(Nmesh, Ndeg);SQPSPYErr1 = zeros(Nmesh, Ndeg);
PUPXErrInf = zeros(Nmesh, Ndeg);PUPSErrInf = zeros(Nmesh, Ndeg);PVPSErrInf = zeros(Nmesh, Ndeg);
PUPXErr2 = zeros(Nmesh, Ndeg);PUPSErr2 = zeros(Nmesh, Ndeg);PVPSErr2 = zeros(Nmesh, Ndeg);
PUPXErr1 = zeros(Nmesh, Ndeg);PUPSErr1 = zeros(Nmesh, Ndeg);PVPSErr1 = zeros(Nmesh, Ndeg);
PVPYErrInf = zeros(Nmesh, Ndeg);PWPSErrInf = zeros(Nmesh, Ndeg);
PVPYErr2 = zeros(Nmesh, Ndeg);PWPSErr2 = zeros(Nmesh, Ndeg);
PVPYErr1 = zeros(Nmesh, Ndeg);PWPSErr1 = zeros(Nmesh, Ndeg);

for n = 1:Ndeg
    for m = 1:Nmesh
        Solver = NonhydroStaticSolver3dTest(Order(n), Order(n), M(m), Mz(m));
        len(m, n) = 20/M(m);
        dofs(m,n) = numel(Solver.meshUnion(1).x);
        PostProcess = NdgPostProcess(Solver.meshUnion(1),strcat('NonhydroStaticSolver3dTest/3d','/','NonhydroStaticSolver3dTest'));
        ExactValue = cell(1);
        ExactValue{1}(:,:,1) = Solver.ExactPSPX;
        ExactValue{1}(:,:,2) = Solver.ExactPSPY;
        ExactValue{1}(:,:,3) = Solver.ExactSPSPX;
        ExactValue{1}(:,:,4) = Solver.ExactSPSPY;
        ExactValue{1}(:,:,5) = Solver.ExactSQPSPX;
        ExactValue{1}(:,:,6) = Solver.ExactSQPSPY;
        ExactValue{1}(:,:,7) = Solver.ExactPUPX;
        ExactValue{1}(:,:,8) = Solver.ExactPUPS;
        ExactValue{1}(:,:,9) = Solver.ExactPVPY;
        ExactValue{1}(:,:,10) = Solver.ExactPVPS;
        ExactValue{1}(:,:,11) = Solver.ExactPWPS;
        
        fphys = cell(1);
        fphys{1}(:,:,1) = Solver.NonhydrostaticSolver.PSPX;
        fphys{1}(:,:,2) = Solver.NonhydrostaticSolver.PSPY;
        fphys{1}(:,:,3) = Solver.NonhydrostaticSolver.SPSPX;
        fphys{1}(:,:,4) = Solver.NonhydrostaticSolver.SPSPY;
        fphys{1}(:,:,5) = Solver.NonhydrostaticSolver.SQPSPX;
        fphys{1}(:,:,6) = Solver.NonhydrostaticSolver.SQPSPY;
        fphys{1}(:,:,7) = Solver.NonhydrostaticSolver.PUPX;
        fphys{1}(:,:,8) = Solver.NonhydrostaticSolver.PUPS;
        fphys{1}(:,:,9) = Solver.NonhydrostaticSolver.PVPY;
        fphys{1}(:,:,10) = Solver.NonhydrostaticSolver.PVPS;
        fphys{1}(:,:,11) = Solver.NonhydrostaticSolver.PWPS;
        
        err = PostProcess.evaluateNormErrInf( fphys, ExactValue );
        PSPXErrInf( m, n ) = err(1);PSPYErrInf( m, n ) = err(2);SPSPXErrInf( m, n ) = err(3);
        SPSPYErrInf( m, n ) = err(4);SQPSPXErrInf( m, n ) = err(5);SQPSPYErrInf( m, n ) = err(6);
        PUPXErrInf( m, n ) = err(7); PUPSErrInf( m, n ) = err(8); PVPYErrInf( m, n ) = err(9); 
        PVPSErrInf( m, n ) = err(10);PWPSErrInf( m, n ) = err(11);
        
        err = PostProcess.evaluateNormErr2( fphys, ExactValue );
        PSPXErr2( m, n ) = err(1);PSPYErr2( m, n ) = err(2);SPSPXErr2( m, n ) = err(3);
        SPSPYErr2( m, n ) = err(4);SQPSPXErr2( m, n ) = err(5);SQPSPYErr2( m, n ) = err(6);
        PUPXErr2( m, n ) = err(7); PUPSErr2( m, n ) = err(8); PVPYErr2( m, n ) = err(9); 
        PVPSErr2( m, n ) = err(10);PWPSErr2( m, n ) = err(11);
        
        err = PostProcess.evaluateNormErr1( fphys, ExactValue );
        PSPXErr1( m, n ) = err(1);PSPYErr1( m, n ) = err(2);SPSPXErr1( m, n ) = err(3);
        SPSPYErr1( m, n ) = err(4);SQPSPXErr1( m, n ) = err(5);SQPSPYErr1( m, n ) = err(6);
        PUPXErr1( m, n ) = err(7); PUPSErr1( m, n ) = err(8); PVPYErr1( m, n ) = err(9); 
        PVPSErr1( m, n ) = err(10);PWPSErr1( m, n ) = err(11);
        
        clear Solver;
        clear PostProcess;
    end
    % print table
    % print table
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), PSPXErr1(:, n), PSPXErr2(:, n), PSPXErrInf(:, n))
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), PSPYErr1(:, n), PSPYErr2(:, n), PSPYErrInf(:, n))
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), SPSPXErr1(:, n), SPSPXErr2(:, n), SPSPXErrInf(:, n))
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), SPSPYErr1(:, n), SPSPYErr2(:, n), SPSPYErrInf(:, n))    
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), SQPSPXErr1(:, n), SQPSPXErr2(:, n), SQPSPXErrInf(:, n))
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), SQPSPYErr1(:, n), SQPSPYErr2(:, n), SQPSPYErrInf(:, n))
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), PUPXErr1(:, n), PUPXErr2(:, n), PUPXErrInf(:, n))
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), PUPSErr1(:, n), PUPSErr2(:, n), PUPSErrInf(:, n)) 
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), PVPYErr1(:, n), PVPYErr2(:, n), PVPYErrInf(:, n))     

    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), PVPSErr1(:, n), PVPSErr2(:, n), PVPSErrInf(:, n))            
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), PWPSErr1(:, n), PWPSErr2(:, n), PWPSErrInf(:, n))  
    
    % plot figure
    co = color{n};
    figure(1); 
    hold on;
   plot(dofs(:, n), PSPXErr1(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PSPYErr1(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SPSPXErr1(:, n), [co, marker{3}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SPSPYErr1(:, n), [co, marker{4}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SQPSPXErr1(:, n), [co, marker{5}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SQPSPYErr1(:, n), [co, marker{6}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PUPXErr1(:, n), [co, marker{7}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PUPSErr1(:, n), [co, marker{8}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PVPYErr1(:, n), [co, marker{9}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PVPSErr1(:, n), [co, marker{10}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PWPSErr1(:, n), [co, marker{11}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );    
    
    figure(2);
    hold on;
   plot(dofs(:, n), PSPXErr2(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PSPYErr2(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SPSPXErr2(:, n), [co, marker{3}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SPSPYErr2(:, n), [co, marker{4}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SQPSPXErr2(:, n), [co, marker{5}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SQPSPYErr2(:, n), [co, marker{6}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PUPXErr2(:, n), [co, marker{7}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PUPSErr2(:, n), [co, marker{8}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PVPYErr2(:, n), [co, marker{9}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PVPSErr2(:, n), [co, marker{10}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PWPSErr2(:, n), [co, marker{11}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );  
    
    figure(3); 
    hold on;
   plot(dofs(:, n), PSPXErrInf(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PSPYErrInf(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SPSPXErrInf(:, n), [co, marker{3}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SPSPYErrInf(:, n), [co, marker{4}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SQPSPXErrInf(:, n), [co, marker{5}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), SQPSPYErrInf(:, n), [co, marker{6}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PUPXErrInf(:, n), [co, marker{7}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PUPSErrInf(:, n), [co, marker{8}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PVPYErrInf(:, n), [co, marker{9}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PVPSErrInf(:, n), [co, marker{10}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
   plot(dofs(:, n), PWPSErrInf(:, n), [co, marker{11}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );      
    
end

ylabel_str = {'$L_1$', '$L_2$', '$L_\infty$'};
fontsize = 12;
for n = 1:3
    h = figure(n);
    box on; grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
%     lendstr = {'p1','p2'};
%     legend(lendstr,'Interpreter','Latex');
    %     columnlegend(2,lendstr, 12);
    
    xlabel('$DOFs$', 'Interpreter', 'Latex', 'FontSize', fontsize);
    ylabel(ylabel_str{n}, 'Interpreter', 'Latex', 'FontSize', fontsize);
end

% time = time./time(end);
% figure(4);
% hold on;



% figure(4);
% hold on;
% for n = 1:Ndeg
%     co = color{n};
%     plot(time(:, n)./max(max(time)), Err2(:, n), [co, marker{1}, linestyle],...
%         'LineWidth', linewidth, ...
%         'MarkerSize', markersize ...
%         );
% end
%
% set(gca, 'XScale', 'log', 'YScale', 'log');
%
% lendstr = {'$H$','$HU$','$HV$'};
% legend(lendstr,'Interpreter','Latex');
% %     columnlegend(2,lendstr, 12);
%
% xlabel('$time \;\rm {(s)}$', 'Interpreter', 'Latex', 'FontSize', fontsize);
% ylabel('$L_2$', 'Interpreter', 'Latex', 'FontSize', fontsize);
% box on; grid on;
end

function t1 = convergence_table(len, err1, err2, errInf)
t1 = table;
t1.len = len(:);
t1.('err1') = err1(:);
t1.('a1') = get_ratio(len, err1);

t1.('err2') = err2(:);
t1.('a2') = get_ratio(len, err2);

t1.('errf') = errInf(:);
t1.('af') = get_ratio(len, errInf);
end

function a = get_ratio(len, err)
Nmesh = numel(len);

a = zeros(Nmesh, 1);
for m = 2:Nmesh
    scal_ratio = log2( len(m)/len(m-1) );
    a(m) = log2( err(m)/err(m-1) )./scal_ratio;
end
end

function columnlegend( numcolumns, str, fontsize, location)
%   Author: Simon Henin <shenin@gc.cuny.edu>
%   Revised by Bill Shawn <bill_shawn@foxmail.com>
if nargin < 4
    location = 'NorthEast';
end
[legend_h,object_h,~,~] = legend( str );
numlines = length(str);
numpercolumn = ceil(numlines/numcolumns);
pos = get(legend_h, 'position');
width = numcolumns*pos(3);
rescale = pos(3)/width;
xdata = get(object_h(numlines+1), 'xdata');
ydata1 = get(object_h(numlines+1), 'ydata');
ydata2 = get(object_h(numlines+3), 'ydata');
sheight = ydata1(1)-ydata2(1);
height = ydata1(1);
line_width = (xdata(2)-xdata(1))*rescale;
spacer = xdata(1)*rescale;
loci = get(gca, 'position');
set(legend_h, 'position', [loci(1) pos(2) width pos(4)]);
col = -1;
for i=1:numlines
    if numpercolumn>1
        if mod(i,numpercolumn)==1
            col = col+1;
        end
    else
        col=i-1;
    end
    if i==1
        linenum = i+numlines;
    else
        linenum = linenum+2;
    end
    labelnum = i;
    position = mod(i,numpercolumn);
    if position == 0
        position = numpercolumn;
    end
    set(object_h(linenum), 'ydata', ...
        [(height-(position-1)*sheight) (height-(position-1)*sheight)]);
    set(object_h(linenum), 'xdata', ...
        [col/numcolumns+spacer col/numcolumns+spacer+line_width]);
    set(object_h(linenum+1), 'ydata', ...
        [height-(position-1)*sheight height-(position-1)*sheight]);
    set(object_h(linenum+1), 'xdata', ...
        [col/numcolumns+spacer*3.5 col/numcolumns+spacer*3.5]);
    set(object_h(labelnum), 'position', ...
        [col/numcolumns+spacer*2+line_width height-(position-1)*sheight]);
end
set(legend_h, 'Color', 'None', 'Box', 'off');
pos = get(legend_h, 'position');
fig_pos = get(gca, 'position');
switch lower(location)
    case {'northeast'}
        set(legend_h, 'position', [pos(1)+fig_pos(3)-pos(3) pos(2) pos(3) pos(4)]);
    case {'southeast'}
        set(legend_h, 'position', [pos(1)+fig_pos(3)-pos(3) fig_pos(2)-pos(4)/2+pos(4)/4 pos(3) pos(4)]);
    case {'southwest'}
        set(legend_h, 'position', [fig_pos(1) fig_pos(2)-pos(4)/2+pos(4)/4 pos(3) pos(4)]);
end
set(legend_h, 'Interpreter', 'Latex','Fontsize',fontsize);
end

