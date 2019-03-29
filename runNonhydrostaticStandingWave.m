function runNonhydrostaticStandingWave
deltax = [2 1.5 1 0.5];
Order = [1 2];
len = deltax;
type = enumStdCell.Quad;
Nmesh = numel(deltax);
Ndeg = numel(Order);
dofs = zeros(Nmesh, Ndeg);
time = zeros(Nmesh, Ndeg);

for n = 1:Ndeg
    for m = 1:Nmesh
        dofs(m, n) = (30 * 6) / (deltax(m) * deltax(m)) * (Order(n)+1).^2;
    end
end

errInf = zeros(Nmesh, Ndeg);
err2 = zeros(Nmesh, Ndeg);
err1 = zeros(Nmesh, Ndeg);

linewidth = 1.5; 
markersize = 8;
c = 'b';
linestyle = '-';

color = {c, c, c, c};
marker = {'o', 's', '^', '*'};

% Error = zeros(numel(deltax)*numel(Order),3);
% EtaError = zeros(numel(deltax), numel(Order));
% UError = zeros(numel(deltax), numel(Order));
% wsError = zeros(numel(deltax), numel(Order));
for n = 1:Ndeg
    for m = 1:Nmesh
  
    Solver = NonhydrostaticStandingWave2d(Order(n),deltax(m),type);
    tic; Solver.matSolve(); time(m, n) = toc;
    PostProcess = NdgPostProcess(Solver.meshUnion(1),strcat('NonhydrostaticStandingWave2d','/','NonhydrostaticStandingWave2d'));
    fext = cell(1);
    fext{1}(:,:,1) = Solver.fexact;fext{1}(:,:,2) = zeros(size(Solver.fext));fext{1}(:,:,3) = zeros(size(Solver.fext));
    fphys = cell(1);
    fphys{1}(:,:,1) = Solver.fphys{1}(:,:,1) - Solver.d;fphys{1}(:,:,2) = zeros(size(Solver.fext));fphys{1}(:,:,3) = zeros(size(Solver.fext));
    error2 = PostProcess.evaluateNormErr2( fphys, fext );
    err2(m, n) = error2(1);
    
    error1 = PostProcess.evaluateNormErr1( fphys, fext );
    err1(m, n) = error1(1);
    
    errorInf = PostProcess.evaluateNormErrInf( fphys, fext );
    errInf(m, n) = errorInf(1);
 
    clear Solver;
    clear PostProcess;
    end
    % print table
    fprintf('\n==================deg = %d==================\n', n);
    convergence_table(len, err1(:, n), err2(:, n), errInf(:, n), ...
        time(:, n))
    
    % plot figure
    co = color{n}; ma = marker{n};
    figure(1); plot(dofs(:, n), err1(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    figure(2); plot(dofs(:, n), err2(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    figure(3); plot(dofs(:, n), errInf(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;
    
    figure(4); plot(dofs(:, n), time(:, n), [co, ma, linestyle],...
        'LineWidth', linewidth, ...
        'MarkerSize', markersize ...
        ); 
    hold on;    
    
end

ylabel_str = {'$L_1$', '$L_2$', '$L_\infty$', '$time(s)$'};
fontsize = 16;
for n = 1:4
    figure(n);
    box on; grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
    lendstr = cell(Ndeg, 1);
    for m = 1:Ndeg
        lendstr{m} = ['$p=', num2str(m), '$'];
    end
    legend(lendstr, 'box', 'off',...
        'Interpreter', 'Latex', 'FontSize', fontsize);
    xlabel('$DOFs$', 'Interpreter', 'Latex', 'FontSize', fontsize);
    ylabel(ylabel_str{n}, 'Interpreter', 'Latex', 'FontSize', fontsize);
end
end

function t1 = convergence_table(len, err1, err2, errInf, time)
t1 = table;
t1.len = len(:);
t1.('err1') = err1(:);
t1.('a1') = get_ratio(len, err1);

t1.('err2') = err2(:);
t1.('a2') = get_ratio(len, err2);

t1.('errf') = errInf(:);
t1.('af') = get_ratio(len, errInf);

t1.('time') = time(:);
end

function a = get_ratio(len, err)
Nmesh = numel(len);

a = zeros(Nmesh, 1);
for m = 2:Nmesh
    scal_ratio = log2( len(m)/len(m-1) );
    a(m) = log2( err(m)/err(m-1) )./scal_ratio;
end
end