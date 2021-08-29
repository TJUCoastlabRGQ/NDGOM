function runEllipticProblemTest3dNew
M = [1/8 1/16 1/32 1/64];
Mz = [2 4 8 16];
Order = [1];

Nmesh = numel(M);
Ndeg = numel(Order);
len = zeros(Nmesh, Ndeg);
dofs = zeros(Nmesh, Ndeg);

linewidth = 1.5;
markersize = 8;

linestyle = '-';

color = {'k', 'r'};  %black for order one, red for order two
marker = {'o'};%circle for the computed field

ErrInf = zeros(Nmesh, Ndeg);
Err2 = zeros(Nmesh, Ndeg);
Err1 = zeros(Nmesh, Ndeg);
for n = 1:Ndeg
    for m = 1:Nmesh
        Solver = EllipticProblemTest3dNew(Order(n), Order(n), M(m), Mz(m));
        Solver.EllipticProblemSolve;
        len(m, n) = 1/Mz(m);
        dofs(m,n) = numel(Solver.meshUnion(1).x);
        PostProcess = NdgPostProcess(Solver.meshUnion(1),strcat('Result/EllipticProblemTest3dNew/3d','/','EllipticProblemTest3dNew'));
        ExactValue = cell(1);
        ExactValue{1} = Solver.ExactSolution;
        fphys = cell(1);
        fphys{1}(:,:,1) = reshape(Solver.SimulatedSolution, Solver.meshUnion.cell.Np, Solver.meshUnion.K);
        for i = 2:PostProcess.Nvar
            fphys{1}(:,:,i) = zeros(size(fphys{1}(:,:,1)));
            ExactValue{1}(:,:,i) = zeros(size(fphys{1}(:,:,1)));
        end
        
        % For all Newmann boundary
        if (all(Solver.meshUnion.BoundaryEdge.ftype == enumBoundaryCondition.SlipWall) && ...
                strcmp(Solver.SurfaceBoundaryEdgeType, 'Newmann') && strcmp(Solver.BottomBoundaryEdgeType,'Newmann'))
            fphys{1}(:,:,1) = fphys{1}(:,:,1) - (fphys{1}(1) - ExactValue{1}(1));
        end
        
        err = PostProcess.evaluateNormErrInf( fphys, ExactValue );
        ErrInf( m, n ) = err(1);
        
        err = PostProcess.evaluateNormErr2( fphys, ExactValue );
        Err2( m, n ) = err(1);
        
        err = PostProcess.evaluateNormErr1( fphys, ExactValue );
        Err1( m, n ) = err(1);
        
        clear Solver;
        clear PostProcess;
    end
    % print table
    % print table
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len(:,n), Err1(:, n), Err2(:, n), ErrInf(:, n))
    
    % plot figure
    co = color{n};
    figure(1); plot(dofs(:, n), Err1(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    hold on;
    
    
    figure(2); plot(dofs(:, n), Err2(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    hold on;
    
    figure(3); plot(dofs(:, n), ErrInf(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );
    hold on;
end

ylabel_str = {'$L_1$', '$L_2$', '$L_\infty$'};
fontsize = 12;
for n = 1:3
    h = figure(n);
    box on; grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    
    lendstr = {'p1','p2'};
    legend(lendstr,'Interpreter','Latex');
    %     columnlegend(2,lendstr, 12);
    
    xlabel('$DOFs$', 'Interpreter', 'Latex', 'FontSize', fontsize);
    ylabel(ylabel_str{n}, 'Interpreter', 'Latex', 'FontSize', fontsize);
end

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

