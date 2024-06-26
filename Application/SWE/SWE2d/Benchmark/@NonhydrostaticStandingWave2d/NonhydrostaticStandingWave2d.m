classdef NonhydrostaticStandingWave2d < SWEPreBlanaced2d
    %NONHYDROSTATICSTANDINGWAVE2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        dt
%         d = 7.612
%         d = 50
        d = 10
        fexact
        A = 0.1;
        T
        Lambda = 20
    end
    
    methods
        function obj = NonhydrostaticStandingWave2d(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEPreBlanaced2d();
            obj.hmin = 1e-3;      
            obj.initPhysFromOptions( mesh );
                             
            c = sqrt(obj.gra*obj.Lambda/2/pi*tanh(2*pi*obj.d/obj.Lambda));
            obj.T = obj.Lambda/c;
            obj.fexact = obj.A * cos(2*pi/obj.Lambda*mesh.x)*cos(2*pi/obj.T*20);
        end
                
        function NonhydroPostprocess(obj)  
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat('Result/',mfilename,'/2d/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Eta = zeros( Ntime,1 );
            exactEta = zeros( Ntime,1 );
            x0 = 18.0;
            h = obj.d;
            a = obj.A;
%             c = sqrt( obj.gra*obj.Lambda/2/pi*tanh(2*pi*h/obj.Lambda) );
%             Period = obj.Lambda / c;
            for t = 1:Ntime
%                 exactEta(t) = -obj.A * cos(2*pi*x0/obj.Lambda + 2*pi*outputTime(t) /obj.T);
                exactEta(t) = obj.A * cos(2*pi/obj.Lambda*x0)*cos(2*pi/obj.T*outputTime(t));
                tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  x0, 0.5, x0, t )-h;
                Eta(t) = tempdata(1);
            end
            figure;
            set(gcf,'position',[50,50,1050,400]);
            plot(outputTime,Eta,'k','LineWidth',1.5);
            hold on;
            set(gca,'YLIM',[-1.3*a, 1.3*a],'Fontsize',15);
            xlabel({'$t\;\rm{(s)}$'},'Interpreter','latex');
            ylabel({'$\eta\;\rm{(m)}$'},'Interpreter','latex');
            
%             str = strcat('Hydro',num2str(obj.d),'.fig');
%             h = openfig(str,'reuse'); % open figure
%             D1=get(gca,'Children'); %get the handle of the line object
%             XData1=get(D1,'XData'); %get the x data
%             YData1=get(D1,'YData'); %get the y data
%             close(h);
%             plot(XData1, YData1,'k--','LineWidth',1.5);
            plot(outputTime(1:10:end),exactEta(1:10:end),'ro','markersize',4.5,'Linewidth',1.5);
            h = legend('SWE','Theory');
            set(h,'fontsize',15);
            legend('boxoff');
        end
        
        function HydroPostprocess(obj)  
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/',mfilename));
            Ntime = PostProcess.Nt;
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
            Eta = zeros( Ntime,1 );
            x0 = 10;
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
%                 fphys{m}(:,:,1) =  obj.A * cos(2*pi*mesh.x/obj.Lambda) * cos( 2*pi * sqrt(obj.gra*obj.d)/obj.Lambda*0) - fphys{m}(:,:,4);
%                 fphys{m}(:,:,1) =  -obj.A * cos( 2*pi*mesh.x/obj.Lambda ) - fphys{m}(:,:,4);  
                fphys{m}(:,:,1) =  obj.A * cos( 2*pi*mesh.x/obj.Lambda ) - fphys{m}(:,:,4);  
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;                
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Hydrostatic;
            option('outputNcfileNum') = 1;
            option('outputType') = enumOutputFile.NetCDF;
%             option('nonhydrostaticType') = enumNonhydrostaticType.Hydrostatic;
        end
    end
    
end

function [ mesh ] = makeUniformMesh(N, deltax, type)
bctype = [...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall];
% bctype = [...
%     enumBoundaryCondition.Clamped, ...
%     enumBoundaryCondition.Clamped, ...
%     enumBoundaryCondition.Clamped, ...
%     enumBoundaryCondition.Clamped];

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 20], [0, 1], 20/deltax, 10, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N,[0, 20], [0, 2*deltax], 20/deltax, 2, bctype);
%     mesh = makeUniformQuadMesh(N,[ -1, 1 ], [ -1, 1 ], deltax, deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
