classdef BohaiEstuary < SWEBarotropic3d
    %BOHAIESTUARY 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties( Constant )
        
        hcrit = 0.1;
        
        SMSFile = [pwd,'\Application\SWE\SWE3d\Benchmark\@BohaiEstuary\Data\fort.14'];
        
        tidalFile = [pwd,'\Application\SWE\SWE3d\Benchmark\@BohaiEstuary\Data\OutputBoundary.xlsx'];
        
        GotmFile = [pwd,'\Application\SWE\SWE3d\Benchmark\@BohaiEstuary\Data\gotmturb.nml'];
    end
    
    properties
        finalTime = 1363*600 + 54*3600;
        
        tideinterval = 600;
        
        OBEdgeIndex
        
        Tide
    end
    
    properties
        T0
        
        S0
        
        alphaT
        
        betaS
    end
    
    methods
        function obj = BohaiEstuary(N, Nz, Mz)
            %BOHAIESTUARY 构造此类的实例
            %   此处显示详细说明
            [ mesh2d, mesh3d ] = obj.makeChannelMesh( N, Nz, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 4];
            obj.outputFieldOrder3d = [ 1 2 4 6 ];
            obj.T0 = 0;
            obj.S0 = 0;
            obj.alphaT = 0;
            obj.betaS = 0;
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.ReadTideElevation;
            %             obj.Tide{1}(:,:,1436) = 0;
        end
        
        drawTopographyAndGaugePoints( obj );
        
        function BohaiEstuaryPostProcess( obj )
            Time = ncread('D:\Sharewithpc\研究工作\20220912\Result\2d\BohaiEstuary.1-1.1.nc','time');
%             N0=datenum([0 12 7 20 10 0]);
            N0=datenum([0 12 8 8 10 0]);
            Newtime = zeros(numel(Time),1);
            Newtime(1) = N0;
            for i = 2:numel(Time)
                Newtime(i) = Newtime(i-1) + datenum([0 0 0 0 0 Time(i) - Time(i-1)]);
            end     
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') );
            EtaIndex = [1 2 3 4 5];
            for i = 1:numel(EtaIndex)
                file = dir([path,'/Data/P',num2str(i),'/*.csv']);
                Eta = zeros(numel(file),1);
                for j = 1:numel(file)
                    Data = xlsread([path,'/Data/P',num2str(i),'/P',num2str(i),'.',num2str(j-1),'.csv']);
                    Eta(j) = Data(6);
                end
                figure;
                hold on;
                Index = Newtime >= datenum([0 12 8 18 0 0]);
                plot(Newtime(Index), Eta(Index),'k','linewidth',1.5);
                MData = xlsread([path,'/Data/Measured-Data/','潮位',num2str(i),'.xlsx']);
                N1 = datenum([0 12 8 21 0 0]);
                MeasureTime = zeros(numel(MData(:,1)),1);
                MeasureTime(1) = N1 + datenum([0 0 0 MData(1,1)*6 0 0]);
                for j = 2:numel(MeasureTime)
                    MeasureTime(j) = MeasureTime(j-1) + datenum([0 0 0 0 0 (MData(j,1) - MData(j-1,1))*6*3600]);
                end
                plot(MeasureTime, MData(:,2),'ro','markersize',7.5,'Linewidth',1.5);
                legend({'Simulation','Measurement'},'Fontsize',12);
                legend('boxoff');
                box on;
                ylabel('$\eta\;(m)$','interpreter','Latex','Fontsize',14);                
                set(gca,'XTickLabelRotation',30);
                datetick('x','mm/dd HH:MM');
                set(gca,'Fontsize',12,'Linewidth',1.4);
                set(gcf,'position',[50,50,700,350]);
            end
            
%             UIndex = [3 4 5];
%             for i = 1:numel(UIndex)
%                 file = dir([path,'/Data/P',num2str(UIndex(i)),'/*.csv']);
%                 Vel = zeros(numel(file),1);
%                 Dir = zeros(numel(file),1);
%                 for j = 1:numel(file)
%                     Data = xlsread([path,'/Data/P',num2str(UIndex(i)),'/P',num2str(UIndex(i)),'.',num2str(j-1),'.csv']);
%                     Vel(j) = Data(1);
%                     if ( Data(7)<0 && Data(8)>=0)
%                         Dir(j) = ...
%                             450 - atan2d(Data(8), Data(7));
%                     else
%                         Dir(j) = ...
%                             90 - atan2d(Data(8), Data(7));
%                     end
%                 end
%                 figure;
%                 hold on;
%                 Index = Newtime >= datenum([0 12 8 18 0 0]);
%                 plot(Newtime(Index), Vel(Index),'k','linewidth',1.5);                            
%                 MData = xlsread([path,'/Data/Measured-Data/','流速',num2str(UIndex(i)),'.xlsx']);
%                 N1 = datenum([0 12 8 21 0 0]);
%                 MeasureTime = zeros(numel(MData(:,1)),1);
%                 MeasureTime(1) = N1 + datenum([0 0 0 MData(1,1)*6 0 0]);
%                 for j = 2:numel(MeasureTime)
%                     MeasureTime(j) = MeasureTime(j-1) + datenum([0 0 0 0 0 (MData(j,1) - MData(j-1,1))*6*3600]);
%                 end
%                 plot(MeasureTime, MData(:,2),'ro','markersize',7.5,'Linewidth',1.5);
%                 legend({'Simulation','Measurement'},'Fontsize',12);
%                 legend('boxoff');
%                 box on;
%                 ylim([0,1.0]);
%                 ylabel('$Velocity\;(m/s)$','interpreter','Latex','Fontsize',14);
%                 set(gca,'XTickLabelRotation',30);
%                 datetick('x','mm/dd HH:MM');
%                 set(gca,'Fontsize',12,'Linewidth',1.4);
%                 set(gcf,'position',[50,50,700,350]);
%                 
%                 figure;
%                 hold on;
%                 plot(Newtime(Index), Dir(Index),'k','linewidth',1.5); 
%                 MData = xlsread([path,'/Data/Measured-Data/','流向',num2str(UIndex(i)),'.xlsx']);
%                 N1 = datenum([0 12 8 21 0 0]);
%                 MeasureTime = zeros(numel(MData(:,1)),1);
%                 MeasureTime(1) = N1 + datenum([0 0 0 MData(1,1)*6 0 0]);
%                 for j = 2:numel(MeasureTime)
%                     MeasureTime(j) = MeasureTime(j-1) + datenum([0 0 0 0 0 (MData(j,1) - MData(j-1,1))*6*3600]);
%                 end  
%                 plot(MeasureTime, MData(:,2),'ro','markersize',7.5,'Linewidth',1.5);
%                 legend({'Simulation','Measurement'},'Fontsize',12);
%                 legend('boxoff');
%                 box on;
%                 ylabel('$Deg$','interpreter','Latex','Fontsize',14);
%                 set(gca,'XTickLabelRotation',30);
%                 datetick('x','mm/dd HH:MM');
%                 set(gca,'Fontsize',12,'Linewidth',1.4);
%                 set(gcf,'position',[50,50,700,350]);
%             end
            
        end
        
    end
    
    methods( Access = protected)
        
        ReadTideElevation( obj );
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
            delta = obj.tideinterval;
            for m1 = 1:obj.Nmesh
                %     mesh = obj.meshUnion( m1 );
                
                % time step
                s1 = floor( time/delta ) + 1;
                s2 = s1 + 1;
                alpha1 = ( delta * (s2 - 1) - time ) / delta;
                alpha2 = ( time - delta * (s1 - 1) ) / delta;
                
                ind = obj.OBEdgeIndex{m1};
                fnT = obj.Tide{m1}(:, :, s1) .* alpha1 + obj.Tide{m1}(:, :, s2) .* alpha2;
                obj.fext2d{m1}(:,ind,1) = max( fnT - obj.fext2d{m1}(:,ind,4), 0 );
                obj.fext3d{m1}(:,:,3) = obj.meshUnion.BoundaryEdge.matVerticalRepmatFacialValue(obj.fext2d{m1}(:,:,1));
            end
        end
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            fphys2d{1} = zeros( obj.mesh2d(1).cell.Np, obj.mesh2d(1).K, obj.Nfield2d );
            fphys2d = obj.matInterpolateTopography(fphys2d);
            fphys2d{1}(:,:,1) = -1*fphys2d{1}(:,:,4);
            
            %             Index = find(obj.meshUnion.mesh2d.x >500005 & obj.meshUnion.mesh2d.y < 4306120);
            %             fphys2d{1}(3*obj.meshUnion.mesh2d.K*obj.meshUnion.mesh2d.cell.Np + Index) = -28.0;
            
            
            %             Index = find(abs(obj.meshUnion.mesh2d.x - 498357)<1 & abs(obj.meshUnion.mesh2d.y - 4311231)<1);
            %             fphys2d{1}(3*obj.meshUnion.mesh2d.K*obj.meshUnion.mesh2d.cell.Np + Index) = -18.0;
            
            fphys{1} = zeros( obj.meshUnion(1).cell.Np, obj.meshUnion(1).K, obj.Nfield );
        end
        
        function fphys2d = matInterpolateTopography( obj, fphys2d )
            [ path, ~, ~ ] = fileparts( mfilename('fullpath') );
            for m = 1:obj.Nmesh
                mesh = obj.mesh2d(m);
                tpfile = [path, '/Data/Bathmetry.txt'];
                fid1 = fopen(tpfile, 'r');
                % READ FIRST LINE
                fgetl(fid1);
                % READ POINT
                temp = fscanf(fid1,'%d',2); Nv = temp(2);
                data = fscanf(fid1,'%d %f %f %f\n',[4,Nv]);
                fclose(fid1);
                interp = scatteredInterpolant( ...
                    data(2,:)', data(3,:)', -1*data(4,:)', 'linear');
                
                fphys2d{m}(:,:,4) = interp(mesh.x, mesh.y);
            end
        end% func
        
        function [ option ] = setOption( obj, option )
            ftime = obj.finalTime;
            startTime = 1363*600;
            outputIntervalNum = 324;
            option('startTime') = startTime;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = (ftime - startTime)/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 230;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('limiterType') = enumLimiter.Vert;
            option('EosType') = enumEOSType.Linear;
            option('ConstantVerticalEddyViscosityValue') = 0.1;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.Constant;
            option('ConstantHorizontalEddyViscosityValue') = 70;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.Beta;
            option('f0 for beta coriolis solver') = 2*sin((38 + 54/60)/180*pi)*7.29*10^(-5);
            %             option('f0 for beta coriolis solver') = 0;
            option('beta for beta coriolis solver') = 0.0;
            option('PhysicalSurfaceRoughnessLength') = 0.0001;
            option('PhysicalBottomRoughnessLength') = 0.0001;
        end
        
        
        function  [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, Mz)
            [ mesh2d ] = makeSMSFileUMeshUnion2d( N, obj.SMSFile );
            cell = StdPrismTri( N, Nz );
            zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
            mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
            mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
            mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
            mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
            mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
            mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
        end
    end
end

