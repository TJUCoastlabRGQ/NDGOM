function runNonhydrostaticSolitaryWave1d
% deltax = [0.25 0.2 0.125 0.1 0.05 0.025];
% deltax = [ 1 0.625 0.5 0.4 0.25 0.2 0.125 0.1 0.05 0.025];
deltax = [ 0.5 0.25 0.125 0.0625];
% deltax = [ 1 0.5];
% deltax = 0.025;
Order =[ 1 2 ];
% Order = 2;
len = deltax;
Nmesh = numel(deltax);
Ndeg = numel(Order);
dofs = zeros(Nmesh, Ndeg);

linewidth = 1.5; 
markersize = 8;

linestyle = '-';

color = {'k', 'r'};  %black for order one, red for order two
marker = {'o', 's', '^'};%circle for Eta, square for U, triangle for Ws
time = zeros(Nmesh, Ndeg);

for n = 1:Ndeg
    for m = 1:Nmesh
        dofs(m, n) = ceil( (50 / deltax(m)) ) * (Order(n)+1);
    end
end

EtaErrInf = zeros(Nmesh, Ndeg); UErrInf = zeros(Nmesh, Ndeg);  WsErrInf = zeros(Nmesh, Ndeg);PErrInf = zeros(Nmesh, Ndeg);
EtaErr2 = zeros(Nmesh, Ndeg);   UErr2 = zeros(Nmesh, Ndeg);    WsErr2 = zeros(Nmesh, Ndeg);PErr2 = zeros(Nmesh, Ndeg);
EtaErr1 = zeros(Nmesh, Ndeg);   UErr1 = zeros(Nmesh, Ndeg);    WsErr1 = zeros(Nmesh, Ndeg);PErr1 = zeros(Nmesh, Ndeg);
for n = 1:Ndeg
    for m = 1:Nmesh
    Solver = NonhydrostaticSolitaryWave1d(Order(n),deltax(m));
    tic;
    Solver.matSolve;
    time(m,n) = toc;
    PostProcess = NdgPostProcess(Solver.meshUnion(1),strcat('NonhydrostaticSolitaryWave1d','/','NonhydrostaticSolitaryWave1d'));
    fext = cell(1);
    fext{1}(:,:,1) = Solver.H5;fext{1}(:,:,2) = Solver.U5;fext{1}(:,:,3) = Solver.W5;fext{1}(:,:,4) = Solver.P5;fext{1}(:,:,5) = zeros( size( Solver.P5 ) );
    fphys = cell(1);

    
    fphys{1}(:,:,1) = Solver.fphys{1}(:,:,1);fphys{1}(:,:,2) = Solver.fphys{1}(:,:,2)./Solver.fphys{1}(:,:,1);
    fphys{1}(:,:,3) = Solver.fphys{1}(:,:,5)./Solver.fphys{1}(:,:,1);fphys{1}(:,:,4) = Solver.fphys{1}(:,:,6);
    fphys{1}(:,:,5) = zeros( size( fphys{1}(:,:,3) ) );
    err = PostProcess.evaluateNormErrInf( fphys, fext );
    EtaErrInf( m, n ) = err(1); UErrInf( m, n ) = err(2); WsErrInf( m, n ) = err(3);PErrInf( m, n ) = err(4);

    err = PostProcess.evaluateNormErr2( fphys, fext );
    EtaErr2( m, n ) = err(1); UErr2( m, n ) = err(2); WsErr2( m, n ) = err(3);PErr2( m, n ) = err(4);    
    
    err = PostProcess.evaluateNormErr1( fphys, fext );
    EtaErr1( m, n ) = err(1); UErr1( m, n ) = err(2); WsErr1( m, n ) = err(3);PErr1( m, n ) = err(4);   
    
    clear Solver;
    clear PostProcess;
    end
    
    % print table
    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len, EtaErr1(:, n), EtaErr2(:, n), EtaErrInf(:, n))

    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len, UErr1(:, n), UErr2(:, n), UErrInf(:, n))   

    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len, WsErr1(:, n), WsErr2(:, n), WsErrInf(:, n)) 

    fprintf('\n==================deg = %d==================\n', Order(n));
    convergence_table(len, PErr1(:, n), PErr2(:, n), PErrInf(:, n)) 
    % plot figure
    co = color{n};
    figure(1); plot(dofs(:, n), EtaErr1(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    plot(dofs(:, n), UErr1(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    
    plot(dofs(:, n), WsErr1(:, n), [co, marker{3}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    
    
    figure(2); plot(dofs(:, n), EtaErr2(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    plot(dofs(:, n), UErr2(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    
    plot(dofs(:, n), WsErr2(:, n), [co, marker{3}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );  
    
    figure(3); plot(dofs(:, n), EtaErrInf(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;    
    plot(dofs(:, n), UErrInf(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    
    plot(dofs(:, n), WsErrInf(:, n), [co, marker{3}, linestyle],...
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

    lendstr = {'$\eta p1$','$U p1$','$W p1$',...
        '$\eta p2$','$U p2$','$W p2$'};
    
    columnlegend(2,lendstr, 12);

    xlabel('$DOFs$', 'Interpreter', 'Latex', 'FontSize', fontsize);
    ylabel(ylabel_str{n}, 'Interpreter', 'Latex', 'FontSize', fontsize);
end

% time = time./time(end);
% figure(4);
% hold on;



    figure(4);
    hold on;
    for n = 1:Ndeg
    co = color{n};
    plot(time(:, n)./max(max(time)), EtaErr2(:, n), [co, marker{1}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    plot(time(:, n)./max(max(time)), UErr2(:, n), [co, marker{2}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    
    plot(time(:, n)./max(max(time)), WsErr2(:, n), [co, marker{3}, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        );  
    end

    set(gca, 'XScale', 'log', 'YScale', 'log');

    lendstr = {'$\eta p1$','$U p1$','$W p1$',...
        '$\eta p2$','$U p2$','$W p2$'};
    
    columnlegend(2,lendstr, 12);
   
 xlabel('$time \;\rm {(s)}$', 'Interpreter', 'Latex', 'FontSize', fontsize);
 ylabel('$L_2$', 'Interpreter', 'Latex', 'FontSize', fontsize);
    box on; grid on;
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
for i=1:numlines,
    if numpercolumn>1
        if mod(i,numpercolumn)==1,
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
    if position == 0,
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
switch lower(location),
    case {'northeast'}
        set(legend_h, 'position', [pos(1)+fig_pos(3)-pos(3) pos(2) pos(3) pos(4)]);
    case {'southeast'}
        set(legend_h, 'position', [pos(1)+fig_pos(3)-pos(3) fig_pos(2)-pos(4)/2+pos(4)/4 pos(3) pos(4)]);
    case {'southwest'}
        set(legend_h, 'position', [fig_pos(1) fig_pos(2)-pos(4)/2+pos(4)/4 pos(3) pos(4)]);
end
set(legend_h, 'Interpreter', 'Latex','Fontsize',fontsize);
end