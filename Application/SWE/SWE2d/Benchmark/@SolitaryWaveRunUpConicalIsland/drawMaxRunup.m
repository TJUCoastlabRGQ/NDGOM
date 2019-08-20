function drawMaxRunup( obj )

% bctype = [...
%     enumBoundaryCondition.ZeroGrad, ...
%     enumBoundaryCondition.ZeroGrad, ...
%     enumBoundaryCondition.ClampedDepth, ...
%     enumBoundaryCondition.ZeroGrad];
%
% % M = 150;
% N = 1;
% mesh = makeUniformQuadMesh(N, [0, 25], [0, 30], 250, 300, bctype);
% %mesh = obj.meshUnion;
%
% % casename = 'ConicalLand2d_150_1';
% % conpos = NdgPostProcess( mesh, casename );

[ path, name, ext ] = fileparts( mfilename('fullpath') );
pathstr = strcat(path, '\result\case 0.045\Maximum Runup.csv');
MeasuredData = xlsread(pathstr);

figure;
pax = polaraxes;
polarplot(MeasuredData(:,1)./360*2*pi,MeasuredData(:,2),'ro');
pax.ThetaZeroLocation = 'bottom';

hold on; 
conpos = NdgPostProcess(obj.meshUnion(1),strcat('SolitaryWaveRunUpConicalIsland','/','SolitaryWaveRunUpConicalIsland'));
dryFlag = zeros( obj.meshUnion(1).cell.Np, obj.meshUnion(1).K );
% threshold = 2.3e-3;
for t = 1:conpos.Nt
    [ fphys ] = conpos.accessOutputResultAtStepNum( t );
    dep = fphys{1}(:,:,1);
%     dryFlag = dryFlag + (dep > 2.5 * obj.hmin); %case 0.096
%     dryFlag = dryFlag + (dep > 2 * obj.hmin); %case 0.181
    dryFlag = dryFlag + (dep > 2.8 * obj.hmin); %case 0.181
end

Xcenter = 12.96;
Ycenter = 13.80;

Index = ( dryFlag == 0 );
Xcor = obj.meshUnion(1).x(Index);
Ycor = obj.meshUnion(1).y(Index);
UniqueY = unique(Ycor);
data = zeros(2*numel(UniqueY), 2);
num = 1;
for i = 1:numel(UniqueY)
    flag = ( Ycor == UniqueY(i) );
    NewX = sort( Xcor(flag), 'ascend');
    %only take the fist and last point
    NewX(2:end-1) = [];
    for j = 1:numel(NewX)
        %distance
        data(num,2) = sqrt( (NewX(j) - Xcenter)^2 + (UniqueY(i) - Ycenter)^2 );
        %angle calculated according to Law of Cosines
        squareA = ( NewX(j) - Xcenter - 0 )^2 + ( UniqueY(i) - Ycenter + 1 )^2;
        squareB = ( NewX(j) - Xcenter )^2 + ( UniqueY(i) - Ycenter )^2;
        squareC = 1;
        
        if NewX(j) >= Xcenter
           data(num,1)  = acos( ( squareB + squareC - squareA ) /(2*sqrt(squareB)*1));
        else
           data(num,1) = 2*pi - acos( ( squareB + squareC - squareA ) /(2*sqrt(squareB)*1)); 
        end
        
        num = num + 1;
    end
    
end

data = sortrows(data);

% index = [1 2 6 7 8 9 12 14 16 18 22 29 32 34 37 38 39 40 42 43 44 45 46 ...
%     48 52 53 60 65 67 69 73 77 79 81 84 88 89 90]; %case 0.096

% index = [1 2 4 10 16 21 22 23 26 33 35 36 38 39 41 43 44 47 50 55 59 72 79 82 86]; %case 0.181

index = [1 4 6 13 17 20 21 22 23 24 25 26 27 28 29 30 33 38 42 43 46 48 50 53 55 58 61 63 66 70 72 74 78 81 86 87 91 93 94];
polarplot(data(index,1),(data(index,2) - 1.1)/1.22,'k','Linewidth',1.5);

set(gca,'Fontsize',12);
% title('$H/H_0=0.096$', 'Interpreter', 'latex','FontSize', 15)
% title('$H/H_0=0.181$', 'Interpreter', 'latex','FontSize', 15)
title('$H/H_0=0.045$', 'Interpreter', 'latex','FontSize', 15)
% legend({'Measured max run-up point', 'Computed dry mesh nodes with $p=1$'}, ...
%     'Interpreter', 'latex', 'FontSize', 16, 'box', 'off', ...
%     'Location', 'northoutside');
% xlabel('$x$ (m)', 'Interpreter', 'latex','FontSize', 16);
% ylabel('$y$ (m)', 'Interpreter', 'latex','FontSize', 16);
end

