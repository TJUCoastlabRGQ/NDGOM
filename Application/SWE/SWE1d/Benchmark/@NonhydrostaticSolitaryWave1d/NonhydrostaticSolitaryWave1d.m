classdef NonhydrostaticSolitaryWave1d < SWEConventional1d
    %NONHYDROSTATICSOLITARYWAVE �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    properties(Constant)
        %> wet/dry depth threshold
        hmin = 0.01
        %> gravity acceleration
        gra = 9.8        
    end    
    
    properties
        H0
        H5
        U0
        U5
        W0
        W5
        P0
        P5
    end
    
    methods
        function obj = NonhydrostaticSolitaryWave1d(N, deltax)
            [ mesh ] = makeUniformMesh( N, deltax);
            obj = obj@SWEConventional1d();
            obj.Solitarywave(mesh);  
            obj.initPhysFromOptions( mesh );
%             obj.matSolve;
%             obj.Postprocess;
        end
        
       function VisualPostprocess(obj)
           Visual = makeVisualizationFromNdgPhys( obj );
           Visual.drawResult( obj.fphys{1}(:, :, 4) );
           Visual.drawHandle.FaceColor = [0.5 0.5 0.5];
           Visual.drawResult( obj.fphys{1}(:, :, 1) );
           Visual.drawHandle.FaceAlpha = 0.95;
           light('Position',[-1 0 0],'Style','local');
        end
        
        function Postprocess(obj)
            mesh = obj.meshUnion(1);
            deltapoint = (mesh.cell.N + 1);
            d = 1;
            a = 0.2;          
            time = [4 8 12 16];
            Nintp = 400;
            xd = linspace(-10, 70, Nintp)';
            yd = 0.025*ones(size(xd));
%             PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/',mfilename));
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
%% plot eta
            figure;
            hold on;            
            set(gcf,'position',[50,50,1050,400]);           
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H0(1:deltapoint:numel(mesh.x)) - d,'k','Linewidth',1.5);            
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H0(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            for i = 1:numel(time);
                [~,Index] = sort(abs(outputTime-time(i)));
                [ fg ] = PostProcess.interpolateOutputStepResultToGaugePoint(  xd, yd, xd, Index(1) );
                plot(xd,fg(:,1)' - d,'k','Linewidth',1.5);
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H4(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H8(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H12(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.H16(1:deltapoint:numel(mesh.x)) - d,'ro','markersize',4.5);
            set(gca,'Fontsize',12);
            xlabel({'$\it x \;\rm{(m)}$'},'Interpreter','latex');
            ylabel({'$\eta \;\rm {(m)}$'},'Interpreter','latex');    
            box on;
%             set(gca,'YLim',[-0.05 a+0.05],'Fontsize',15);   
%% plot U
            figure;
            hold on;            
            set(gcf,'position',[50,50,1050,400]);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U0(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U0(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            for i = 1:numel(time);           
                [~,Index] = sort(abs(outputTime-time(i)));
                [ fg ] = PostProcess.interpolateOutputStepResultToGaugePoint(  xd, yd, xd, Index(1) );
                plot(xd,fg(:,2)'./fg(:,1)','k','Linewidth',1.5);
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U4(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U8(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U12(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U16(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            set(gca,'Fontsize',12);
            xlabel({'$\it x(m)$'},'Interpreter','latex');
            ylabel({'$U\;\rm{(m/s)}$'},'Interpreter','latex');    
            box on;
%             set(gca,'YLim',[-0.05 a],'Fontsize',15);                   
%% plot W
            figure;
            hold on;            
            set(gcf,'position',[50,50,1050,400]);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W0(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W0(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            for i = 1:numel(time);                      
                [~,Index] = sort(abs(outputTime-time(i)));
                [ fg ] = PostProcess.interpolateOutputStepResultToGaugePoint(  xd, yd, xd, Index(1) );
                plot(xd,fg(:,4)'./fg(:,1)','k','Linewidth',1.5);
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W4(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W8(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W12(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W16(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            set(gca,'Fontsize',12);
            xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
            ylabel({'$w\;\rm {(m/s)}$'},'Interpreter','latex');   
            box on;
%             set(gca,'YLim',[-0.3*a 0.3*a+0.01],'Fontsize',15);       
%% plot P
            figure;
            hold on;            
            set(gcf,'position',[50,50,1050,400]);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.P0(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.P0(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            for i = 1:numel(time);                      
                [~,Index] = sort(abs(outputTime-time(i)));
                [ fg ] = PostProcess.interpolateOutputStepResultToGaugePoint(  xd, yd, xd, Index(1) );
                plot(xd,fg(:,5)','k','Linewidth',1.5);
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.P4(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.P8(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.P12(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.P16(1:deltapoint:numel(mesh.x)),'ro','markersize',4.5);
            set(gca,'Fontsize',12);
            xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
            ylabel({'$P$'},'Interpreter','latex');   
            box on;
        end
    end
    
    methods(Access = protected)
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,1) = obj.H0;
%                 fphys{m}(:,:,1) = obj.Eta0;
                fphys{m}(:,:,2) = fphys{m}(:,:,1).*obj.U0;
                fphys{m}(:,:,5) =  fphys{m}(:,:,1).*obj.W0;
                fphys{m}(:,:,6) =  obj.P0;
%                 fphys{m}(:,:,6) =  fphys{m}(:,:,1).*obj.W0./2;
            end       
        end
        
        
        function [ option ] = setOption( obj, option )
            ftime = 5;
            outputIntervalNum = 4500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('outputFieldOrder') = [1 2 3 6 7];
            option('limiterType') = enumLimiter.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
%             option('nonhydrostaticType') = enumNonhydrostaticType.Hydrostatic;
        end
    end
    
end

function [ mesh ] = makeUniformMesh( N, deltax)
xlim = [-10, 30];
M = ceil(( xlim(2) - xlim(1) ) / deltax);
bcType = [...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
mesh.InnerEdge.Ne = mesh.InnerEdge.Ne + 1;
mesh.BoundaryEdge.ftype = [];
mesh.BoundaryEdge.Ne = [];
FToV = zeros(1, numel(mesh.InnerEdge.FToV) + 1 );
FToV(1) = mesh.BoundaryEdge.FToV(1);
% FToV(end) = mesh.BoundaryEdge.FToV(2);
FToV(2:end) = mesh.InnerEdge.FToV;
mesh.InnerEdge.FToV = FToV;
mesh.BoundaryEdge.FToV = [];

FToE = zeros(2, mesh.InnerEdge.Ne);
FToE(1) = mesh.K; FToE(2) = 1;
FToE(3:end) = mesh.InnerEdge.FToE(:);
mesh.InnerEdge.FToE = FToE;
mesh.BoundaryEdge.FToE = [];

FToF = zeros(2, mesh.InnerEdge.Ne);
FToF(1) = 2; FToF(2) = 1;
FToF(3:end) = mesh.InnerEdge.FToF(:);
mesh.InnerEdge.FToF = FToF;
mesh.BoundaryEdge.FToF = [];

FToN1 = zeros(1,mesh.InnerEdge.Ne);
FToN2 = zeros(1,mesh.InnerEdge.Ne);
FToN1(1) = mesh.BoundaryEdge.FToN1(2);
FToN1(2:end) = mesh.InnerEdge.FToN1;
mesh.InnerEdge.FToN1 = FToN1;
FToN2(1) = mesh.BoundaryEdge.FToN1(1);
FToN2(2:end) = mesh.InnerEdge.FToN2;
mesh.InnerEdge.FToN2 = FToN2;
mesh.BoundaryEdge.FToN1 = [];
mesh.BoundaryEdge.FToN2 = [];

mesh.InnerEdge.nx = ones(1, mesh.InnerEdge.Ne);
mesh.BoundaryEdge.nx = [];
mesh.InnerEdge.ny = zeros(1, mesh.InnerEdge.Ne);
mesh.BoundaryEdge.ny = [];
mesh.InnerEdge.nz = zeros(1, mesh.InnerEdge.Ne);
mesh.BoundaryEdge.nz = [];
mesh.InnerEdge.Js = ones(1, mesh.InnerEdge.Ne);
mesh.BoundaryEdge.Js = [];

end% func