obj = load('Solver.mat');
obj = obj.Solver;
PostProcess = NdgPostProcess(obj.meshUnion(1),'NonhydrostaticStandingWave2d/NonhydrostaticStandingWave2d');
Ntime = PostProcess.Nt;
outputTime = ncread( PostProcess.outputFile{1}, 'time' );
Entropy = zeros( Ntime,1 ); ExactEntropy = zeros( Ntime,1 );
Energy = zeros( Ntime,1 ); ExactEnergy = zeros( Ntime,1 );
syms z;

for t = 1:Ntime
    fphys = PostProcess.accessOutputResultAtStepNum(t);
    %     averageUA = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * (fphys{1}(:,:,2)./fphys{1}(:,:,1))) );
    %     averageVA = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * (fphys{1}(:,:,3)./fphys{1}(:,:,1))) );
    %     averageHA = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * fphys{1}(:,:,1)) );
    %     Entropy(t) = sum( obj.gra * averageHA ) + sum( (averageUA.^2 + averageVA.^2)./obj.meshUnion(1).LAV./2 );
    %     Energy(t) = sum( 0.5 * obj.gra * averageHA.^2./obj.meshUnion(1).LAV + 0.5 * averageHA .* ( averageUA.^2 + averageVA.^2 )./(obj.meshUnion(1).LAV.^2) );
    
    %     [ ExactEta, ExactU ] = GetExactSolution( obj, outputTime(t) );
    %     averageExactUA = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * ExactU) );
    %     averageExactHA = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * (ExactEta + obj.d)) );
    %     ExactEntropy(t) = sum( obj.gra * averageHA ) + sum( averageExactUA.^2 ./obj.meshUnion(1).LAV./2 );
    %     ExactEnergy(t) = sum( 0.5 * obj.gra * averageExactHA.^2./obj.meshUnion(1).LAV + 0.5 * averageExactHA .* ( averageExactUA.^2  )./(obj.meshUnion(1).LAV.^2) );
    
    
    averageU2A = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * (fphys{1}(:,:,2).^2./fphys{1}(:,:,1).^2)) );
    averageV2A = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * (fphys{1}(:,:,3).^2./fphys{1}(:,:,1).^2)) );
    averageHA = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * fphys{1}(:,:,1)) );
    averageH2A = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * fphys{1}(:,:,1).^2) );
    averageHU2V2A = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * ( fphys{1}(:,:,2).^2./fphys{1}(:,:,1) + fphys{1}(:,:,3).^2./fphys{1}(:,:,1)) ) );
    Entropy(t) = sum( obj.gra * averageHA ) + sum( (averageU2A + averageV2A)./2 );
    Energy(t) = sum( 0.5 * obj.gra * averageH2A + 0.5 * averageHU2V2A );
    
    %     [ ExactEta, ExactU ] = GetExactSolution( obj, outputTime(t) );
    ExactEta = obj.A * cos(pi*obj.meshUnion(1).x/obj.Lambda) * cos( pi * sqrt(obj.gra*obj.d)/obj.Lambda*outputTime(t));
    ExactU = obj.A * sqrt(obj.gra * obj.d)/obj.d*sin(pi*obj.meshUnion(1).x/obj.Lambda) *  sin( pi * sqrt(obj.gra*obj.d)/obj.Lambda*outputTime(t));
%     A = obj.A; T = obj.T; Lambda = obj.Lambda; d = obj.d;
%     x = obj.meshUnion(1).x;
%     parfor i =  1:numel(ExactU)
%         display(i);
%         ExactU(i) = int ( A * 2 * pi /T * cosh( 2*pi/Lambda * (z + d) ) ./ sinh (2*pi/Lambda*d) .*...
%             sin( 2*pi/Lambda.*x(i))*sin( 2*pi/T*outputTime(t) ),z, -d , ExactEta(i) )./(ExactEta(i) + d);
%     end
    
    averageExactHA = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * ( ExactEta + obj.d )) );
    averageExactU2A = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * ( ExactU.^2) ) );
    averageExactH2A = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * ( (ExactEta + obj.d).^2 )) );
    averageExactHU2A = sum( obj.meshUnion(1).J .* ( obj.meshUnion(1).cell.M * ( (ExactEta + obj.d) .* ExactU.^2) ) );
    
    ExactEntropy(t) = sum( obj.gra * averageExactHA ) + sum( averageExactU2A./2 );
    ExactEnergy(t) = sum( 0.5 * obj.gra * averageExactH2A + 0.5 * averageExactHU2A );
    
    %     tempEntropy = obj.gra * fphys{1}(:,:,1) + (fphys{1}(:,:,2).^2 + fphys{1}(:,:,3).^2)./(fphys{1}(:,:,1).^2)./2;
    %     tempEnergy = 0.5 * obj.gra .* fphys{1}(:,:,1).^2 + fphys{1}(:,:,1) .* (fphys{1}(:,:,2).^2 + fphys{1}(:,:,3).^2)./(fphys{1}(:,:,1).^2)./2;
    %     for k = 1:obj.meshUnion(1).K
    %         tempEntropy(:,k) =  obj.meshUnion(1).cell.V \ tempEntropy(:,k);
    %         tempEnergy(:,k) =  obj.meshUnion(1).cell.V \ tempEnergy(:,k);
    %     end
    %     Entropy(t) = sum( tempEntropy(1,:) .* obj.meshUnion(1).LAV );
    %     Energy(t) = sum( tempEnergy(1,:) .* obj.meshUnion(1).LAV );
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

% function [ ExactEta, ExactU ] = GetExactSolution( obj, outputTime )
% syms z;
% mesh = obj.meshUnion( 1 );
% ExactEta = obj.A * cos( 2*pi/obj.Lambda*mesh.x ) .* cos( 2*pi/obj.T*outputTime );
% ExactU = zeros(size(ExactEta));
% parfor i =  1:numel(ExactU)
% ExactU(i) = int ( obj.A * 2 * pi /obj.T * cosh( 2*pi/obj.Lambda * (z + obj.d) ) ./ sinh (2*pi/obj.Lambda*obj.d) .* sin( 2*pi/obj.Lambda.*mesh.x(i))*sin( 2*pi/obj.T*outputTime ),z, -obj.d , ExactEta(i) )./(ExactEta(i) + obj.d);
% end
% end