classdef NonhydrostaticStandingWave1d < SWEConventional1d
    %NONHYDROSTATICSTANDINGWAVE2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        rho = 1000
        %> wet/dry depth threshold
        hmin = 0.01
        %> gravity acceleration
        gra = 9.8        
    end
    
    properties
        dt
        d = 7.612
%         d = 50
        fexact
        A = 0.1;
        T
        Lambda = 20
    end
    
    methods
        function obj = NonhydrostaticStandingWave1d(N, M)
            [ mesh ] = makeUniformMesh( N, M );
            obj = obj@SWEConventional1d();
%             obj.outputFieldOrder = [1 2 3 6];   
            obj.initPhysFromOptions( mesh );
                             
            c = sqrt(obj.gra*obj.Lambda/2/pi*tanh(2*pi*obj.d/obj.Lambda));
            obj.T = obj.Lambda/c;
            obj.fexact = -obj.A * cos(2*pi/obj.Lambda*mesh.x)*cos(2*pi/obj.T*10);
        end
                
        function NonhydroPostprocess(obj)  
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Eta = zeros( Ntime,1 );
            exactEta = zeros( Ntime,1 );
            x0 = 12.8;
            h = obj.d;
            a = obj.A;
%             c = sqrt(obj.gra*obj.Lambda/2/pi*tanh(2*pi*h/obj.Lambda));
            for t = 1:Ntime
                exactEta(t) = - obj.A * cos(2*pi/obj.Lambda*x0)*cos(2*pi/obj.T*outputTime(t));
                tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  x0, 0, 0, t )-h;
                Eta(t) = tempdata(1);
            end
            figure;
            set(gcf,'position',[50,50,1050,400]);
            plot(outputTime,Eta,'k','LineWidth',1.5);
            hold on;
            set(gca,'YLIM',[-1.3*a, 1.3*a],'Fontsize',12);
            xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex');
            ylabel({'$\eta\;\rm{(m)}$'},'Interpreter','latex');
            
%             str = strcat('Hydro',num2str(obj.d),'.fig');
%             h = openfig(str,'reuse'); % open figure
%             D1=get(gca,'Children'); %get the handle of the line object
%             XData1=get(D1,'XData'); %get the x data
%             YData1=get(D1,'YData'); %get the y data
%             close(h);
%             plot(XData1, YData1,'k--','LineWidth',1.5);
            plot(outputTime,exactEta,'ro','markersize',1.5);
            legend('Nonhydro','Hydro','Exact');
            legend('boxoff');
        end
        
        function HydroPostprocess(obj)  
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Eta = zeros( Ntime,1 );
            x0 = 17.5;
            h = obj.d;
            a = obj.A;
            for t = 1:Ntime
                tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  x0, 0.2, x0, t )-h;
                Eta(t) = tempdata(1);
            end
            figure;
            set(gcf,'position',[50,50,1050,400]);
            plot(outputTime,Eta,'k','LineWidth',1.5);
            legend('Hydro');
            legend('boxoff');
            set(gca,'YLIM',[-1*a, 1.5*a],'Fontsize',12);
            xlabel({'$\it t(s)$'},'Interpreter','latex');
            ylabel({'$\eta(m)$'},'Interpreter','latex');
            str = strcat('Hydro',num2str(obj.d),'.fig');
            saveas(gca,str);
            close(gcf);
        end
    end
    
    methods(Access = protected)
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                bot = -obj.d;
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,4) = bot;
%                 fphys{m}(:,:,1) =  obj.d;
                fphys{m}(:,:,1) =  -obj.A * cos(2*pi*mesh.x/obj.Lambda) - fphys{m}(:,:,4);                
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 500;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;                
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('limiterType') = enumLimiter.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
%             option('nonhydrostaticType') = enumNonhydrostaticType.Hydrostatic;
        end
    end
    
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 20];
bcType = [enumBoundaryCondition.SlipWall, enumBoundaryCondition.SlipWall];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
mesh.InnerEdge.Ne = mesh.InnerEdge.Ne + 1;
mesh.BoundaryEdge.ftype = zeros(2,0);
mesh.BoundaryEdge.Ne = 0;
FToV = zeros(1, numel(mesh.InnerEdge.FToV) + 1 );
FToV(1) = mesh.BoundaryEdge.FToV(1);
% FToV(end) = mesh.BoundaryEdge.FToV(2);
FToV(2:end) = mesh.InnerEdge.FToV;
mesh.InnerEdge.FToV = FToV;
mesh.BoundaryEdge.FToV = zeros(0,2);

FToE = zeros(2, mesh.InnerEdge.Ne);
FToE(1) = mesh.K; FToE(2) = 1;
LAV = zeros(1, mesh.InnerEdge.Ne);
LAV(1) = mesh.BoundaryEdge.LAV(1);
LAV(2:end) = mesh.InnerEdge.LAV;
FToE(3:end) = mesh.InnerEdge.FToE(:);
mesh.InnerEdge.FToE = FToE;
mesh.InnerEdge.LAV = LAV;
mesh.BoundaryEdge.FToE = zeros(2,0);
mesh.BoundaryEdge.LAV = zeros(2,0);

FToF = zeros(2, mesh.InnerEdge.Ne);
FToF(1) = 2; FToF(2) = 1;
FToF(3:end) = mesh.InnerEdge.FToF(:);
mesh.InnerEdge.FToF = FToF;
mesh.BoundaryEdge.FToF = zeros(2,0);

FToN1 = zeros(1,mesh.InnerEdge.Ne);
FToN2 = zeros(1,mesh.InnerEdge.Ne);
FToN1(1) = mesh.BoundaryEdge.FToN1(2);
FToN1(2:end) = mesh.InnerEdge.FToN1;
mesh.InnerEdge.FToN1 = FToN1;
FToN2(1) = mesh.BoundaryEdge.FToN1(1);
FToN2(2:end) = mesh.InnerEdge.FToN2;
mesh.InnerEdge.FToN2 = FToN2;
mesh.BoundaryEdge.FToN1 = zeros(2,0);
mesh.BoundaryEdge.FToN2 = zeros(2,0);

mesh.InnerEdge.nx = ones(1, mesh.InnerEdge.Ne);
mesh.BoundaryEdge.nx = zeros(2,0);
mesh.InnerEdge.ny = zeros(1, mesh.InnerEdge.Ne);
mesh.BoundaryEdge.ny = zeros(0,2);
mesh.InnerEdge.nz = zeros(1, mesh.InnerEdge.Ne);
mesh.BoundaryEdge.nz = zeros(0,2);
mesh.InnerEdge.Js = ones(1, mesh.InnerEdge.Ne);
mesh.BoundaryEdge.Js = zeros(0,2);

mesh.EToE(1) = mesh.K;
mesh.EToE(end) = 1;

mesh.EToF(1) = 2;
mesh.EToF(end) = 1;
end% func
