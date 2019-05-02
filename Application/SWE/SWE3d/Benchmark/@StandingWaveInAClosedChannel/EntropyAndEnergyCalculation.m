function EntropyAndEnergyCalculation(obj)
PostProcess = NdgPostProcess(obj.mesh2d(1),'StandingWaveInAClosedChannel/StandingWaveInAClosedChannel');
Ntime = PostProcess.Nt;
outputTime = ncread( PostProcess.outputFile{1}, 'time' );
Entropy = zeros( Ntime,1 );
Energy = zeros( Ntime,1 );

for t = 1:Ntime
    fphys = PostProcess.accessOutputResultAtStepNum(t);
    tempEntropy = obj.gra * fphys{1}(:,:,1) + (fphys{1}(:,:,2).^2 + fphys{1}(:,:,3).^2)./(fphys{1}(:,:,1).^2)./2;
    tempEnergy = 0.5 * obj.gra .* fphys{1}(:,:,1).^2 + fphys{1}(:,:,1) .* (fphys{1}(:,:,2).^2 + fphys{1}(:,:,3).^2)./(fphys{1}(:,:,1).^2)./2;
    for k = 1:obj.mesh2d(1).K
        tempEntropy(:,k) =  obj.mesh2d(1).cell.V \ tempEntropy(:,k);
        tempEnergy(:,k) =  obj.mesh2d(1).cell.V \ tempEnergy(:,k);
    end
    Entropy(t) = sum( tempEntropy(1,:) .* obj.mesh2d(1).LAV );
    Energy(t) = sum( tempEnergy(1,:) .* obj.mesh2d(1).LAV );
end
figure;
set(gcf,'position',[50,50,1050,400]);
plot(outputTime,Entropy,'k','LineWidth',1.5);
hold on;
xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex');
ylabel('Entropy');
legend('Entropy');

figure;
set(gcf,'position',[50,50,1050,400]);
plot(outputTime,Energy,'k','LineWidth',1.5);
hold on;
xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex');
ylabel('Energy');
legend('Energy');
end