obj = load('Solver.mat');
obj = obj.Solver;
PostProcess = NdgPostProcess(obj.mesh2d(1),'StandingWaveInAClosedChannel/StandingWaveInAClosedChannel');
Ntime = PostProcess.Nt;
outputTime = ncread( PostProcess.outputFile{1}, 'time' );
Entropy = zeros( Ntime,1 ); ExactEntropy = zeros( Ntime,1 );
Energy = zeros( Ntime,1 ); ExactEnergy = zeros( Ntime,1 );

for t = 1:Ntime
    fphys = PostProcess.accessOutputResultAtStepNum(t);
%     tempEntropy = obj.gra * fphys{1}(:,:,1) + (fphys{1}(:,:,2).^2 + fphys{1}(:,:,3).^2)./(fphys{1}(:,:,1).^2)./2;
%     tempEnergy = 0.5 * obj.gra .* fphys{1}(:,:,1).^2 + fphys{1}(:,:,1) .* (fphys{1}(:,:,2).^2 + fphys{1}(:,:,3).^2)./(fphys{1}(:,:,1).^2)./2;
%     for k = 1:obj.mesh2d(1).K
%         tempEntropy(:,k) =  obj.mesh2d(1).cell.V \ tempEntropy(:,k);
%         tempEnergy(:,k) =  obj.mesh2d(1).cell.V \ tempEnergy(:,k);
%     end
%     Entropy(t) = sum( tempEntropy(1,:) .* obj.mesh2d(1).LAV );
%     Energy(t) = sum( tempEnergy(1,:) .* obj.mesh2d(1).LAV );

    averageU2A = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * (fphys{1}(:,:,2).^2./fphys{1}(:,:,1).^2)) );
    averageV2A = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * (fphys{1}(:,:,3).^2./fphys{1}(:,:,1).^2)) );
    averageHA = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * fphys{1}(:,:,1)) );
    averageH2A = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * fphys{1}(:,:,1).^2) );
    averageHU2V2A = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * ( fphys{1}(:,:,2).^2./fphys{1}(:,:,1) + fphys{1}(:,:,3).^2./fphys{1}(:,:,1)) ) );
    Entropy(t) = sum( obj.gra * averageHA ) + sum( (averageU2A + averageV2A)./2 );
    Energy(t) = sum( 0.5 * obj.gra * averageH2A + 0.5 * averageHU2V2A );
    
    ExactEta = obj.A * cos(pi*obj.mesh2d(1).x/obj.Lambda) * cos( pi * sqrt(obj.gra*obj.H0)/obj.Lambda*outputTime(t));
    ExactU = obj.A * sqrt(obj.gra * obj.H0)/obj.H0*sin(pi*obj.mesh2d(1).x/obj.Lambda) *  sin( pi * sqrt(obj.gra*obj.H0)/obj.Lambda*outputTime(t));
    
    averageExactHA = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * ( ExactEta + obj.H0 )) );
    averageExactU2A = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * ( ExactU.^2) ) );
    averageExactH2A = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * ( (ExactEta + obj.H0).^2 )) );
    averageExactHU2A = sum( obj.mesh2d(1).J .* ( obj.mesh2d(1).cell.M * ( (ExactEta + obj.H0) .* ExactU.^2) ) );
    
    ExactEntropy(t) = sum( obj.gra * averageExactHA ) + sum( averageExactU2A./2 );
    ExactEnergy(t) = sum( 0.5 * obj.gra * averageExactH2A + 0.5 * averageExactHU2A );    

end
figure;
set(gcf,'position',[50,50,1050,400]);
plot(outputTime,Entropy,'k','LineWidth',1.5);
hold on;
plot(outputTime,ExactEntropy,'r','LineWidth',1.5);
xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex');
ylabel('Entropy');
legend('Simulated','Exact');

figure;
set(gcf,'position',[50,50,1050,400]);
plot(outputTime,Energy,'k','LineWidth',1.5);
hold on;
plot(outputTime,ExactEnergy,'r','LineWidth',1.5);
xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex');
ylabel('Energy');
legend('Simulated','Exact');