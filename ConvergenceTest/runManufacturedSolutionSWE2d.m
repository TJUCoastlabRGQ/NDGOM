function runManufacturedSolutionSWE2d
% M = [ 20 10 5 2.5 ];




M = [10 5 3 2.5 ];
Order = [1 2];

As = 900 * M(1) * 2;



% Order =[ 1 2 3 ];
% len = deltax;
Nmesh = numel(M);
Ndeg = numel(Order);
len = zeros(Nmesh, Ndeg);
dofs = zeros(Nmesh, Ndeg);

linewidth = 1.5;
markersize = 8;

linestyle = '-';

color = {'k', 'r'};  %black for order one, red for order two
marker = {'o','s','^'};%circle for H, square for HU, triangle for HV
time = zeros(Nmesh, Ndeg);

% for n = 1:Ndeg
%     for m = 1:Nmesh
%         dofs(m, n) = ceil( (50 / deltax(m)) ) * (Order(n)+1);
%     end
% end

HErrInf = zeros(Nmesh, Ndeg); HUErrInf = zeros(Nmesh, Ndeg);  HVErrInf = zeros(Nmesh, Ndeg);
HErr2 = zeros(Nmesh, Ndeg);   HUErr2 = zeros(Nmesh, Ndeg);    HVErr2 = zeros(Nmesh, Ndeg);
HErr1 = zeros(Nmesh, Ndeg);   HUErr1 = zeros(Nmesh, Ndeg);    HVErr1 = zeros(Nmesh, Ndeg);
for n = 1:Ndeg
    for m = 1:Nmesh
        Solver = ManufacturedSolution2d(Order(n), M(m), enumStdCell.Quad);
        tic;
        Solver.matSolve;
        A = sum( Solver.meshUnion.LAV );
        time(m,n) = toc;
        len(m, n) = M(m);
        dofs(m,n) = numel(Solver.fphys{1}(:,:,1));
        PostProcess = NdgPostProcess(Solver.meshUnion(1),strcat('ManufacturedSolution2d/2d','/','ManufacturedSolution2d'));
        ExactValue = cell(1);
        ExactValue{1} = Solver.ExactValue{1};
        fphys = cell(1);
        fphys{1}(:,:,1:3) = Solver.fphys{1}(:,:,1:3);
        err = PostProcess.evaluateNormErrInf( fphys, ExactValue );
        
        HErrInf( m, n ) = err(1); HUErrInf( m, n ) = err(2); HVErrInf( m, n ) = err(3);
        
        err = PostProcess.evaluateNormErr2( fphys, ExactValue );
        err = err * sqrt(A/As);
        HErr2( m, n ) = err(1); HUErr2( m, n ) = err(2); HVErr2( m, n ) = err(3);
        
        err = PostProcess.evaluateNormErr1( fphys, ExactValue );
        HErr1( m, n ) = err(1); HUErr1( m, n ) = err(2); HVErr1( m, n ) = err(3);
        
        clear Solver;
        clear PostProcess;
    end
    % print table
    % print table
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), HErr1(:, n), HErr2(:, n), HErrInf(:, n))
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), HUErr1(:, n), HUErr2(:, n), HUErrInf(:, n))
    
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), HVErr1(:, n), HVErr2(:, n), HVErrInf(:, n))
    
    % plot figure
    co = color{n};
    figure(1); plot(dofs(:, n), HErr1(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    hold on;
    plot(dofs(:, n), HUErr1(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    
    plot(dofs(:, n), HVErr1(:, n), [co, marker{3}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    
    
    figure(2); plot(dofs(:, n), HErr2(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    hold on;
    plot(dofs(:, n), HUErr2(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    
    plot(dofs(:, n), HVErr2(:, n), [co, marker{3}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    
    figure(3); plot(dofs(:, n), HErrInf(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    hold on;
    plot(dofs(:, n), HUErrInf(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    
    plot(dofs(:, n), HVErrInf(:, n), [co, marker{3}, linestyle],...
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
    
    lendstr = {'H1',...
        'HU1','HV1','H2',...
        'HU2','HV2'};
%     legend(lendstr,'Interpreter','Latex');
    columnlegend(2,lendstr, 12);
    
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