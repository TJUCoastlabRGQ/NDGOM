function runNonhydrostaticSolitaryWave
% deltax = [0.25 0.2 0.125 0.1 0.05 0.025];
deltax = [ 1 0.625 0.5 0.4 0.25 0.2 0.125 0.1];
Order =[ 1 2 ];
% Order = 2;
len = deltax;
type = enumStdCell.Quad;
Nmesh = numel(deltax);
Ndeg = numel(Order);
dofs = zeros(Nmesh, Ndeg);

linewidth = 1.5; 
markersize = 8;

linestyle = '-';

color = {'k', 'r'};  %black for order one, red for order two
marker = {'o', 's', '^'};%circle for Eta, square for U, triangle for Ws

for n = 1:Ndeg
    for m = 1:Nmesh
        dofs(m, n) = (40 * deltax(m)) / (deltax(m) * deltax(m)) * (Order(n)+1).^2;
    end
end

EtaErrInf = zeros(Nmesh, Ndeg); UErrInf = zeros(Nmesh, Ndeg);  WsErrInf = zeros(Nmesh, Ndeg);PErrInf = zeros(Nmesh, Ndeg);
EtaErr2 = zeros(Nmesh, Ndeg);   UErr2 = zeros(Nmesh, Ndeg);    WsErr2 = zeros(Nmesh, Ndeg);PErr2 = zeros(Nmesh, Ndeg);
EtaErr1 = zeros(Nmesh, Ndeg);   UErr1 = zeros(Nmesh, Ndeg);    WsErr1 = zeros(Nmesh, Ndeg);PErr1 = zeros(Nmesh, Ndeg);
for n = 1:Ndeg
    for m = 1:Nmesh
    Solver = NonhydrostaticSolitaryWave(Order(n),deltax(m),type);
    Solver.matSolve;
    PostProcess = NdgPostProcess(Solver.meshUnion(1),strcat('NonhydrostaticSolitaryWave','/','NonhydrostaticSolitaryWave'));
    fext = cell(1);
    fext{1}(:,:,1) = Solver.H5;fext{1}(:,:,2) = Solver.U5;fext{1}(:,:,3) = Solver.W5;fext{1}(:,:,4) = Solver.P5;fext{1}(:,:,5) = zeros( size( Solver.P5 ) );
    fphys = cell(1);
    %> methods from Lu
%     fphys{1}(:,:,1) = Solver.fphys{1}(:,:,1) - Solver.Depth;fphys{1}(:,:,2) = Solver.fphys{1}(:,:,2)./Solver.fphys{1}(:,:,1);fphys{1}(:,:,3) = Solver.fphys{1}(:,:,6)*2;fphys{1}(:,:,4) = zeros(size(Solver.fphys{1}(:,:,6)));
    %> methods from Stelling
    fphys{1}(:,:,1) = Solver.fphys{1}(:,:,1);fphys{1}(:,:,2) = Solver.fphys{1}(:,:,2)./Solver.fphys{1}(:,:,1);
    fphys{1}(:,:,3) = Solver.fphys{1}(:,:,6)./Solver.fphys{1}(:,:,1);fphys{1}(:,:,4) = Solver.fphys{1}(:,:,7);
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
save EtaErrInf;
save UErrInf;
save WsErrInf;
save PErrInf;
save EtaErr1;
save UErr1;
save WsErr1;
save PErr1;
save EtaErr2;
save UErr2;
save WsErr2;
save PErr2;
ylabel_str = {'$L_1$', '$L_2$', '$L_\infty$'};
fontsize = 16;
for n = 1:3
    figure(n);
    box on; grid on;
    set(gca, 'XScale', 'log', 'YScale', 'log');
%     lendstr = cell(Ndeg, 1);
%     for m = 1:Ndeg
%         lendstr{m} = ['$p=', num2str(m), '$'];
%     end
%     legend(lendstr, 'box', 'off',...
%         'Interpreter', 'Latex', 'FontSize', fontsize);
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