classdef NonhydrostaticSolitaryWave3d < SWEBarotropic3d
    %NONHYDROSTATICSOLITARYWAVE 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        hcrit = 0.01;
        %         ChLength = 100;
        ChLength = 100;
        %> channel width
        ChWidth = 1/3;
        %> channel depth
        InitialDepth = 1;
        %> x range
        %> start time
        startTime = 0;
        %         %> final time
        finalTime = 20;
    end
    
    properties
        H0
        U0
        W0
        CEU0
        CEP0
        CEW0
        H5
        CEU5
        CEP5
        CEW5
        H10
        CEU10
        CEP10
        CEW10
        H15
        CEU15
        CEP15
        CEW15
        H20
        CEU20
        CEP20
        CEW20
    end
    
    properties
        SurfaceBoundaryEdgeType = 'Dirichlet'
    end
    
    
    methods
        function obj = NonhydrostaticSolitaryWave3d(N, Nz, M, Mz)
            [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz);
            obj.Solitarywave(mesh3d);
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [1 2 4 11 13];
            obj.Nfield = 13;
            obj.fieldName3d = {'hu','hv','omega', 'h','nv','z','eta','zx','zy','w', 'hw','hc','p'};
            obj.Nvar = 3;
            obj.varFieldIndex = [1 2 11];
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, mesh3d );
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( mesh2d, mesh3d );
        end
        
        function VisualPostprocess(obj)
            Visual = makeVisualizationFromNdgPhys( obj );
            Visual.drawResult( obj.fphys{1}(:, :, 4) );
            Visual.drawHandle.FaceColor = [0.5 0.5 0.5];
            Visual.drawResult( obj.fphys{1}(:, :, 1) );
            Visual.drawHandle.FaceAlpha = 0.95;
            light('Position',[-1 0 0],'Style','local');
        end
        
        function TimeSerisPostprocess( obj )
            AEtaMax = max(max(obj.H0)) - obj.InitialDepth;
            AUMax = max(max(abs(obj.CEU0)));
            AWMax = max(max(abs(obj.CEW0)));
            APMax = max(max(abs(obj.CEP0)));
            PostProcess = NdgPostProcess(obj.meshUnion(1).mesh2d,strcat('Result/',mfilename,'/2d/',mfilename));
            PostProcess3d = NdgPostProcess(obj.meshUnion(1),strcat('Result/',mfilename,'/3d/',mfilename));
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            etaRatio = zeros(numel(outputTime),1);
            uRatio = zeros(numel(outputTime),1);
            wRatio = zeros(numel(outputTime),1);
            pRatio = zeros(numel(outputTime),1);
            for i = 1:numel(outputTime)
                fresult = PostProcess3d.accessOutputResultAtStepNum( i ); %
                Hu = obj.meshUnion.cell.VCV * fresult{1}(:,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,1);
                Hw = obj.meshUnion.cell.VCV * fresult{1}(:,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,4);
                H = obj.meshUnion.cell.VCV * fresult{1}(:,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,3);
                p = obj.meshUnion.cell.VCV * fresult{1}( :,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,5 );
                uMax = max(max(Hu./H));
                wMax = max(max(Hw./H));
                pMax = max(max(p));
                etaRatio(i) = ( max(max(H)) - obj.InitialDepth - AEtaMax )/0.2 * 100;
                uRatio(i) = ( uMax - AUMax )/ AUMax * 100;
                wRatio(i) = ( wMax - AWMax )/ AWMax * 100;
                pRatio(i) = ( pMax - APMax) /APMax * 100;
            end
            figure;
            set(gcf,'position',[50,50,1050,400]);
            set(gca,'Fontsize',14);
            hold on;
            plot(outputTime,etaRatio,'r','Linewidth',1.5);
            plot(outputTime,uRatio,'g','Linewidth',1.5);
            plot(outputTime,wRatio,'b','Linewidth',1.5);
            plot(outputTime,pRatio,'k','Linewidth',1.5);
            xlabel({'$\it t \;\rm{(s)}$'},'Interpreter','latex');
            ylabel({'RE (%)'},'Interpreter','latex');
%             ylabel({'$ RE $''(%)'},'Interpreter','latex');
            lendstr = {'$\eta$','$u$','$w$',...
                '$p$'};
            box on;
            columnlegend(2,lendstr, 14);
            set(gca,'LineWidth',1.5);
        end
        
        function Postprocess(obj)
            AEtaMax = max(max(obj.H0-obj.InitialDepth));
            AUMax = max(max(abs(obj.CEU0)));
            AWMax = max(max(abs(obj.CEW0)));
            APMax = max(max(abs(obj.CEP0)));
            mesh = obj.meshUnion(1).mesh2d;
            deltapoint = (mesh.cell.N + 1);
            d = 1;
            a = 0.2;
            time = [5 10 15 20];
            Nintp = 400;
            xd = linspace(0, obj.ChLength, Nintp)';
            yd = 1/6*ones(size(xd));
            %             PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            PostProcess = NdgPostProcess(obj.meshUnion(1).mesh2d,strcat('Result/',mfilename,'/2d/',mfilename));
            PostProcess3d = NdgPostProcess(obj.meshUnion(1),strcat('Result/',mfilename,'/3d/',mfilename));
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            %% plot eta
            figure;
            hold on;
            set(gcf,'position',[50,50,1050,400]);
%             plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H0(1:deltapoint:numel(mesh.x)) - d,'k','Linewidth',1.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H0(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            for i = 1:numel(time)
                [~,Index] = sort(abs(outputTime-time(i)));
                [ fg ] = PostProcess.interpolateOutputStepResultToGaugePoint(  xd, yd, xd, Index(1) );
                plot(xd,fg(:,1)' - d,'k','Linewidth',1.5);
                if i == numel(time)
                    EtaMax = max(fg(:,1)) - d;
                    DiffEtaRatio = ( EtaMax - AEtaMax)/AEtaMax * 100;
                    disp(abs(DiffEtaRatio));
                end
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H5(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H10(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H15(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H20(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            %             plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H30(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            text(15 - 5*1.75,0.175,'$0\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*5 - 5*1.75,0.175,'$5\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*10 - 5*1.75,0.175,'$10\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*15 - 5*1.75,0.175,'$15\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*20 - 5*1.75,0.175,'$20\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            
            set(gca,'Fontsize',14);
            xlabel({'$\it x \;\rm{(m)}$'},'Interpreter','latex');
            ylabel({'$\eta \;\rm {(m)}$'},'Interpreter','latex');
            box on;
            set(gca,'LineWidth',1.5);
            
            
            % %             set(gca,'YLim',[-0.05 a+0.05],'Fontsize',15);
            %% plot U
            figure;
            hold on;
            set(gcf,'position',[50,50,1050,400]);
%             plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEU0(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEU0(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            for i = 1:numel(time)
                [~,Index] = sort(abs(outputTime-time(i)));
                fresult = PostProcess3d.accessOutputResultAtStepNum( Index(1) ); %
                Hu = obj.meshUnion.cell.VCV * fresult{1}(:,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,1);
                H = obj.meshUnion.cell.VCV * fresult{1}(:,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,3);
                u = Hu./H;
                plot(mesh.x(1:deltapoint:numel(mesh.x)),u(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
                if i == numel(time)
                    UMax = max(max(u));
                    DiffURatio = (UMax - AUMax)/AUMax * 100;
                    disp(abs(DiffURatio));
                end
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEU5(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEU10(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEU15(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEU20(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            
            text(15 - 5*1.75,0.7*8/10-0.1,'$0\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*5 - 5*1.75,0.7*8/10-0.1,'$5\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*10 - 5*1.75,0.7*8/10-0.1,'$10\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*15 - 5*1.75,0.7*8/10-0.1,'$15\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*20 - 5*1.75,0.7*8/10-0.1,'$20\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            %             plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U16(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            set(gca,'Fontsize',14);
            xlabel({'$\it x(m)$'},'Interpreter','latex');
            ylabel({'$u\;\rm{(m/s)}$'},'Interpreter','latex');
            box on;
            set(gca,'LineWidth',1.5);
            %             set(gca,'YLim',[-0.05 a],'Fontsize',15);
            %% plot W
            figure;
            hold on;
            set(gcf,'position',[50,50,1050,400]);
%             plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEW0(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEW0(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            for i = 1:numel(time)
                [~,Index] = sort(abs(outputTime-time(i)));
                fresult = PostProcess3d.accessOutputResultAtStepNum( Index(1) ); %
                Hv = obj.meshUnion.cell.VCV * fresult{1}(:,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,4);
                H = obj.meshUnion.cell.VCV * fresult{1}(:,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,3);
                v = Hv./H;
                plot(mesh.x(1:deltapoint:numel(mesh.x)),v(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
                if i == numel(time)
                    WMax = max(max(abs(v)));
                    DiffWRatio = (WMax - AWMax)/AWMax * 100;
                    disp(abs(DiffWRatio));
                end
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEW5(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEW10(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEW15(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEW20(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            
            text(15 - 5*1.75,0.2*8/10-0.1,'$0\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*5 - 5*1.75,0.2*8/10-0.1,'$5\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*10 - 5*1.75,0.2*8/10-0.1,'$10\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*15 - 5*1.75,0.2*8/10-0.1,'$15\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*20 - 5*1.75,0.2*8/10-0.1,'$20\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            % %             plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W16(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            set(gca,'Fontsize',14);
            xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
            ylabel({'$w\;\rm {(m/s)}$'},'Interpreter','latex');
            box on;
            set(gca,'LineWidth',1.5);
            % %             set(gca,'YLim',[-0.3*a 0.3*a+0.01],'Fontsize',15);
            %% plot P
            figure;
            hold on;
            set(gcf,'position',[50,50,1050,400]);
%             plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEP0(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEP0(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            for i = 1:numel(time)
                [~,Index] = sort(abs(outputTime-time(i)));
                fresult = PostProcess3d.accessOutputResultAtStepNum( Index(1) ); %
                p = obj.meshUnion.cell.VCV * fresult{1}( :,(1 + obj.meshUnion(1).Nz )/2:obj.meshUnion(1).Nz:end,5 );
                plot(mesh.x(1:deltapoint:numel(mesh.x)),p(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
                if i == numel(time)
                    PMax = max(max(abs(p)));
                    DiffPRatio = (PMax - APMax)/APMax * 100;
                    disp(abs(DiffPRatio));
                end
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEP5(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEP10(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEP15(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.CEP20(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            
            text(15 - 5*1.75,1000*8/10+4800,'$0\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*5 - 5*1.75,1000*8/10+4800,'$5\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*10 - 5*1.75,1000*8/10+4800,'$10\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*15 - 5*1.75,1000*8/10+4800,'$15\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            text(15 + 3.4259*20 - 5*1.75,1000*8/10+4800,'$20\;\rm s\rightarrow$','Interpreter','latex','Fontsize',14);
            % %             plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.P16(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            set(gca,'Fontsize',14);
            xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
            ylabel({'$p\;\rm {(Pa)}$'},'Interpreter','latex');
            set(gca,'LineWidth',1.5);
            box on;
        end
    end
    
    methods(Access = protected)
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                %water depth
                fphys2d{m}(:,:,1) =  obj.H0;
                fphys{m}(:,:,1) =  obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) ) .* obj.U0;
                fphys{m}(:,:,11) =  obj.meshUnion(1).Extend2dField( fphys2d{1}(:, :, 1) ) .* obj.W0;
            end
        end
        
        %         function fphys = setInitialField( obj )
        %             fphys = cell( obj.Nmesh, 1 );
        %             for m = 1:obj.Nmesh
        %                 mesh = obj.meshUnion(m);
        %                 fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
        %                 fphys{m}(:,:,1) = obj.H0 .* obj.U;
        % %                 fphys{m}(:,:,1) = obj.Eta0;
        %                 fphys{m}(:,:,2) = fphys{m}(:,:,1).*obj.U0;
        %                 fphys{m}(:,:,6) =  fphys{m}(:,:,1).*obj.W0;
        %                 fphys{m}(:,:,7) =  obj.P0;
        % %                 fphys{m}(:,:,6) =  fphys{m}(:,:,1).*obj.W0./2;
        %             end
        %         end
        
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 1;
            
        end
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

% mesh2d = makeUniformQuadMesh( N, ...
%     [ 0, obj.ChLength ], [0, obj.ChWidth], M, ceil(obj.ChWidth/(obj.ChLength/M)), bctype);
mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, ceil(obj.ChWidth/(obj.ChLength/M)), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end

function columnlegend( numcolumns, str, fontsize, location)
%   Author: Simon Henin <shenin@gc.cuny.edu>
%   Revised by Bill Shawn <bill_shawn@foxmail.com>
if nargin < 4
    location = 'NorthEast';
end
% [legend_h,object_h,~,~] = legend( str );
[legend_h,object_h,~,~] = legend( str,'Interpreter', 'Latex','Fontsize',fontsize);
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
