classdef NonhydrostaticSolitaryWave < SWEConventional2d
    %NONHYDROSTATICSOLITARYWAVE 此处显示有关此类的摘要
    %   此处显示详细说明
    
  properties(Constant)
        
        rho = 1000
        Depth = 10
    end
    
    properties
        A = 2
        Eta0
        Eta1
        Eta2
        Eta25
        Eta50
        U0
        U1
        U2
        U25
        U50
        W0
        W1
        W2
        W25
        W50
    end
    
    methods
        function obj = NonhydrostaticSolitaryWave(N, deltax, cellType)
            [ mesh ] = makeUniformMesh(N, deltax, cellType);
            obj = obj@SWEConventional2d();
            obj.hmin = 1e-3;
            obj.Solitarywave(mesh);  
            obj.initPhysFromOptions( mesh );       
            obj.matSolve;
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
            deltapoint = 24;
            mesh = obj.meshUnion(1);      
            d = obj.Depth;
            a = obj.A;          
            time = [25 50];
            Nintp = 2000;
            xd = linspace(0, 900, Nintp)';
            yd = 0.25*ones(size(xd));
%             PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'.',num2str(obj.Nmesh),'-','1','/',mfilename));
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/',mfilename));
            outputTime = ncread( PostProcess.outputFile{1}, 'time' );
%% plot eta
            figure;
            hold on;            
            set(gcf,'position',[50,50,1050,400]);           
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.Eta0(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);            
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.Eta0(1:deltapoint:numel(mesh.x)),'ro','markersize',3);
            for i = 1:numel(time);
                [~,Index] = sort(abs(outputTime-time(i)));
                [ fg ] = PostProcess.interpolateOutputStepResultToGaugePoint(  xd, yd, xd, Index(1) );
                plot(xd,fg(:,1)'- d,'k','Linewidth',1.5);
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.Eta25(1:deltapoint:numel(mesh.x)),'ro','markersize',3);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.Eta50(1:deltapoint:numel(mesh.x)),'ro','markersize',3);
            set(gca,'Fontsize',12);
            xlabel({'$\it x \;\rm{(m)}$'},'Interpreter','latex');
            ylabel({'$\eta \;\rm {(m)}$'},'Interpreter','latex');    
            box on;
            set(gca,'YLim',[-0.05 a+0.05],'Fontsize',15);   
%% plot U
            figure;
            hold on;            
            set(gcf,'position',[50,50,1050,400]);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U0(1:deltapoint:numel(mesh.x)),'k','Linewidth',1.5);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U0(1:deltapoint:numel(mesh.x)),'ro','markersize',3);
            for i = 1:numel(time);           
                [~,Index] = sort(abs(outputTime-time(i)));
                [ fg ] = PostProcess.interpolateOutputStepResultToGaugePoint(  xd, yd, xd, Index(1) );
                plot(xd,fg(:,2)'./fg(:,1)','k','Linewidth',1.5);
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U25(1:deltapoint:numel(mesh.x)),'ro','markersize',3);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.U50(1:deltapoint:numel(mesh.x)),'ro','markersize',3);
            set(gca,'Fontsize',12);
            xlabel({'$\it x(m)$'},'Interpreter','latex');
            ylabel({'$U\;\rm{(m/s)}$'},'Interpreter','latex');    
            box on;
            set(gca,'YLim',[-0.05 a],'Fontsize',15);                   
%% plot W
            figure;
            hold on;            
            set(gcf,'position',[50,50,1050,400]);
            plot(mesh.x(1:8:numel(mesh.x)),obj.W0(1:8:numel(mesh.x)),'k','Linewidth',1.5);
            plot(mesh.x(1:8:numel(mesh.x)),obj.W0(1:8:numel(mesh.x)),'ro','markersize',3);
            for i = 1:numel(time);                      
                [~,Index] = sort(abs(outputTime-time(i)));
                [ fg ] = PostProcess.interpolateOutputStepResultToGaugePoint(  xd, yd, xd, Index(1) );
                plot(xd,2*fg(:,4)'./fg(:,1)','k','Linewidth',1.5);
            end
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W25(1:deltapoint:numel(mesh.x)),'ro','markersize',3);
            plot(mesh.x(1:deltapoint:numel(mesh.x)),obj.W50(1:deltapoint:numel(mesh.x)),'ro','markersize',3);            
            set(gca,'Fontsize',12);
            xlabel({'$x\;\rm{(m)}$'},'Interpreter','latex');
            ylabel({'$w_s\;\rm {(m/s)}$'},'Interpreter','latex');   
            box on;
            set(gca,'YLim',[-0.3*a 0.3*a+0.01],'Fontsize',15);       
        end
    end
    
    methods(Access = protected)
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                fphys{m}(:,:,1) = obj.Eta0 + obj.Depth;
%                 fphys{m}(:,:,1) = obj.Eta0;
                fphys{m}(:,:,2) = fphys{m}(:,:,1).*obj.U0;
                fphys{m}(:,:,6) =  obj.W0./2;
            end       
        end
        
        
        function [ option ] = setOption( obj, option )
            ftime = 1;
            outputIntervalNum = 4500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('outputFieldOrder') = [1 2 3 6];
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;
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

if (type == enumStdCell.Tri)
    mesh = makeUniformTriMesh(N, [0, 450], [0, 3], 450/deltax, 3/deltax, bctype);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformQuadMesh(N, [0, 600], [0, 6], 600/deltax, 6/deltax, bctype);
else
    msgID = [mfile, ':inputCellTypeError'];
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func