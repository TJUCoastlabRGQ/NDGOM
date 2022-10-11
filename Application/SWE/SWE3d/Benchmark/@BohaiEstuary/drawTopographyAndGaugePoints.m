function drawTopographyAndGaugePoints( obj )
drawTopography( obj );
drawGaugePoints( obj );
writeTecplotFile( obj );
writeBottomTecplotFile( obj );
end

function drawTopography( obj )
patch(obj.mesh2d.x, obj.mesh2d.y, obj.fphys2d{1}(:,:,4), 'Edgecolor','none');
colormap jet;
box on;
set(gca,'Linewidth',1.5);
set(gca,'Fontsize',14);
xlabel('$x\;(m)$','Interpreter','Latex','Fontsize',14);
ylabel('$y\;(m)$','interpreter','Latex','Fontsize',14);
axis equal;
h = colorbar;
h.Label.String = '$m$';
h.Label.Interpreter = 'Latex';
h.Label.FontSize = 14;
end

function drawGaugePoints( obj )
% close all
patch(obj.mesh2d.x, obj.mesh2d.y, obj.fphys2d{1}(:,:,4), 'Edgecolor','none');
colormap jet;
% patch(obj.mesh2d.x, obj.mesh2d.y, obj.fphys2d{1}(:,:,4),'Facecolor','none','Edgecolor','none');
hold on;
GP = [446553.7409, 4310967.491;
    445601.0536, 4308376.988;
    448744.2658, 4308763.083;
    457687.493, 4313172.592;
    454665.6124, 4303002.943];
% scatter(GP(:,1), GP(:,2), 20, 'filled', 'k');
axis equal;
set(gca,'xlim',[432910 462610],'ylim',[4297500 4320900]);
for i=1:size(GP,1)
    scatter(GP(i), GP(i+5), 20, 'filled', 'k');
    text(GP(i)-1000,GP(i+5)-1000,['G',num2str(i)],'Fontsize',13,'color','k');
end
box on;
set(gca,'Linewidth',1.5);
set(gca,'Fontsize',14);
xlabel('$x\;(m)$','Interpreter','Latex','Fontsize',14);
ylabel('$y\;(m)$','interpreter','Latex','Fontsize',14);
end

function writeTecplotFile( obj )
Time = ncread('D:\Sharewithpc\研究工作\20220912\Result\2d\BohaiEstuary.1-1.1.nc','time');
fphys = ncread('D:\Sharewithpc\研究工作\20220912\Result\3d\BohaiEstuary.1-1.1.nc','fphys');
Np2d = obj.mesh2d.cell.Np;
K2d = obj.mesh2d.K;
x = obj.mesh2d.x;
y = obj.mesh2d.y;
fid1 = fopen('D:\Sharewithpc\研究工作\20221009\Surface.tec','w');
fprintf(fid1,'%s\n','Title = "ProjectPostProcess "');
fprintf(fid1,'%s\n','VARIABLES ="X","Y","Z","h","u","v","zeta"');
fid2 = fopen('D:\Sharewithpc\研究工作\20221009\Middle.tec','w');
fprintf(fid2,'%s\n','Title = "ProjectPostProcess "');
fprintf(fid2,'%s\n','VARIABLES ="X","Y","Z","h","u","v","zeta"');
fid3 = fopen('D:\Sharewithpc\研究工作\20221009\Bottom.tec','w');
fprintf(fid3,'%s\n','Title = "ProjectPostProcess "');
fprintf(fid3,'%s\n','VARIABLES ="X","Y","Z","h","u","v","zeta"');
Tp = [173 214];
for i = 1:numel(Tp)
    Data = fphys(:,:,:,Tp(i));
    %The top layer
    %     hu = obj.meshUnion.cell.VCV * Data(:,1:5:end,1);
    %     hv = obj.meshUnion.cell.VCV * Data(:,1:5:end,2);
    %     h = obj.meshUnion.cell.VCV * Data(:,1:5:end,3);
    %     z = obj.meshUnion.cell.VCV * Data(:,1:5:end,4);
    hu = Data(4:6,1:5:end,1);
    hv = Data(4:6,1:5:end,2);
    h = Data(4:6,1:5:end,3);
    z = Data(4:6,1:5:end,4);
    fprintf(fid1,'%s %s %s %s %s %s\n',['ZONE T ="P_',num2str(Time(Tp(i))),'",F=FEPOINT,ET=TRIANGLE,N=',num2str(Np2d*K2d),' E=',num2str(K2d)]);
    for p = 1:Np2d*K2d
        fprintf(fid1,'%f %f %f %f %f %f %f\n',[x(p) y(p) 0 h(p) hu(p)/h(p) hv(p)/h(p) h(p)+z(p)]);
    end
    for k =1:K2d
        fprintf(fid1,'%d %d %d\n',[(k-1)*3+1 (k-1)*3+2 (k-1)*3+3]);
    end
    % The middle layer
    hu = Data(1:3,3:5:end,1);
    hv = Data(1:3,3:5:end,2);
    h = Data(1:3,3:5:end,3);
    z = Data(1:3,3:5:end,4);
    fprintf(fid2,'%s %s %s %s %s %s\n',['ZONE T ="P_',num2str(Time(Tp(i))),'",F=FEPOINT,ET=TRIANGLE,N=',num2str(Np2d*K2d),' E=',num2str(K2d)]);
    for p = 1:Np2d*K2d
        fprintf(fid2,'%f %f %f %f %f %f %f\n',[x(p) y(p) 0 h(p) hu(p)/h(p) hv(p)/h(p) h(p)+z(p)]);
    end
    for k =1:K2d
        fprintf(fid2,'%d %d %d\n',[(k-1)*3+1 (k-1)*3+2 (k-1)*3+3]);
    end
    % The bottom layer
    %     hu = obj.meshUnion.cell.VCV * Data(:,5:5:end,1);
    %     hv = obj.meshUnion.cell.VCV * Data(:,5:5:end,2);
    %     h = obj.meshUnion.cell.VCV * Data(:,5:5:end,3);
    %     z = obj.meshUnion.cell.VCV * Data(:,5:5:end,4);
    hu = Data(1:3,5:5:end,1);
    hv = Data(1:3,5:5:end,2);
    h = Data(1:3,5:5:end,3);
    z = Data(1:3,5:5:end,4);
    fprintf(fid3,'%s %s %s %s %s %s\n',['ZONE T ="P_',num2str(Time(Tp(i))),'",F=FEPOINT,ET=TRIANGLE,N=',num2str(Np2d*K2d),' E=',num2str(K2d)]);
    for p = 1:Np2d*K2d
        fprintf(fid3,'%f %f %f %f %f %f %f\n',[x(p) y(p) 0 h(p) hu(p)/h(p) hv(p)/h(p) h(p)+z(p)]);
    end
    for k =1:K2d
        fprintf(fid3,'%d %d %d\n',[(k-1)*3+1 (k-1)*3+2 (k-1)*3+3]);
    end
end
fclose(fid1);
fclose(fid2);
fclose(fid3);
end

function writeBottomTecplotFile( obj )
Np2d = obj.mesh2d.cell.Np;
K2d = obj.mesh2d.K;
x = obj.mesh2d.x;
y = obj.mesh2d.y;
fid1 = fopen('D:\Sharewithpc\研究工作\20221009\BottomTopography.tec','w');
fprintf(fid1,'%s\n','Title = "ProjectPostProcess "');
fprintf(fid1,'%s\n','VARIABLES ="X","Y","Z","Bot"');
z = obj.fphys2d{1}(:,:,4);
fprintf(fid1,'%s %s %s %s %s\n',['ZONE T ="P_0','",F=FEPOINT,ET=TRIANGLE,N=',num2str(Np2d*K2d),' E=',num2str(K2d)]);
for p = 1:Np2d*K2d
    fprintf(fid1,'%f %f %f %f\n',[x(p) y(p) 0 z(p)]);
end
for k =1:K2d
    fprintf(fid1,'%d %d %d\n',[(k-1)*3+1 (k-1)*3+2 (k-1)*3+3]);
end
fclose(fid1);
end