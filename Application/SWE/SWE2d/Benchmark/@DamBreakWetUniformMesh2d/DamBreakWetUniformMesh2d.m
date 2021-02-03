classdef DamBreakWetUniformMesh2d < SWEConventional2d
    
    properties( SetAccess = protected )
        theta
    end
    
    properties(Constant)
        %> wet/dry depth threshold
        %         hmin = 1e-4
        %> gravity acceleration
        %         gra = 9.8
        %> Dam position
        damPosition = 500
        %> water depth at the upstream (left)
        h0 = 10
        %> water depth at the downstream (right)
        h1 = 2
    end
    
    methods
        function obj = DamBreakWetUniformMesh2d(N, M, cellType, theta)
            [ mesh ] = makeUniformMesh(N, M, cellType, theta);
            obj = obj@SWEConventional2d();
            obj.hmin = 1e-4;
            obj.theta = theta;
            obj.initPhysFromOptions( mesh );
        end
        
        function CheckSection( obj )
            Ng = 100;
            xg = linspace(0, 1000, Ng)'; yg = zeros(Ng, 1);
            pos = Analysis2d( obj, xg, yg );
            fphyInterp = pos.InterpGaugeResult( obj.fphys );
            ftime = obj.getOption('finalTime');
            %             fext = obj.getExactFunction( ftime );
            fext = getExactFunction(obj, obj.meshUnion(1).cell.Np, obj.meshUnion(1).K,...
                obj.meshUnion(1).x, obj.meshUnion(1).y, ftime);
            fextInterp = pos.InterpGaugeResult( fext );
            
            figure('Color', 'w'); hold on; grid on; box on;
            plot( xg, fphyInterp(:,1), 'b.-' );
            plot( xg, fextInterp(:,1), 'r.-' );
            figure('Color', 'w'); hold on; grid on; box on;
            plot( xg, fphyInterp(:,2), 'b.-' );
            plot( xg, fextInterp(:,2), 'r.-' );
            figure('Color', 'w'); hold on; grid on; box on;
            uphyInterp = fphyInterp(:,2)./fphyInterp(:,1);
            uextInterp = fextInterp(:,2)./fextInterp(:,1);
            plot( xg, uphyInterp, 'b.-' );
            plot( xg, uextInterp, 'r.-' );
        end
        
        function CheckPointTimeHistory( obj )
            xg = 0;
            yg = 0;
            PostProcess = NdgPostProcess(obj.meshUnion(1),strcat(mfilename,'/2d/',mfilename));
            time = PostProcess.time;
            Exacth = zeros(1,PostProcess.Nt);
            Exacthu = zeros(1,PostProcess.Nt);
            Exacthv = zeros(1,PostProcess.Nt);
            Simulatedh = zeros(1,PostProcess.Nt);
            Simulatedhu = zeros(1,PostProcess.Nt);
            Simulatedhv = zeros(1,PostProcess.Nt);
            for i = 1:PostProcess.Nt
                ExactData = getExactFunction(obj, 1, 1,...
                    xg, yg, time{1}(i));
                Exacth(i)  = ExactData{1}(1);
                Exacthu(i) = ExactData{1}(2);
                Exacthv(i) = ExactData{1}(3);
                tempdata = PostProcess.interpolateOutputStepResultToGaugePoint(  xg, yg, 0, i );
                Simulatedh(i) = tempdata(1);
                Simulatedhu(i) = tempdata(2);
                Simulatedhv(i) = tempdata(3);
            end
            figure(1), hold on;
            plot(time{1}, Exacth,'r','Linewidth',1.0);
            plot(time{1}, Simulatedh,'r--','Linewidth',1.0);
            legend({'Theory','Simulation'},'Fontsize',12);
            legend BOXOFF;
            figure(2), hold on;
            plot(time{1}, Exacthu,'r','Linewidth',1.0);
            plot(time{1}, Simulatedhu,'r--','Linewidth',1.0); 
            legend({'Theory','Simulation'},'Fontsize',12);
            legend BOXOFF;
            figure(3), hold on;
            plot(time{1}, Exacthv,'r','Linewidth',1.0);
            plot(time{1}, Simulatedhv,'r--','Linewidth',1.0);   
            legend({'Theory','Simulation'},'Fontsize',12);
            legend BOXOFF;
        end
        
    end
    
    methods(Access=protected)
        function fphys = setInitialField( obj )
            fphys = getExactFunction(obj, obj.meshUnion(1).cell.Np, obj.meshUnion(1).K,...
                obj.meshUnion(1).x, obj.meshUnion(1).y, 0);
        end
        
        function matUpdateExternalField( obj, time, ~ )
            obj.fext = getExactFunction(obj, obj.meshUnion(1).cell.Nfp(1), obj.meshUnion(1).BoundaryEdge.Ne,...
                obj.meshUnion(1).BoundaryEdge.xb, obj.meshUnion(1).BoundaryEdge.yb, time);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 80;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('temporalDiscreteType') = enumTemproalInterval.DeltaTime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            %             option('NumFluxType') = enumSWENumFlux.LF;
        end
        
        fphys = getExactFunction( obj, Np, Ne, x, y, time );
    end
end

function [ mesh, theta ] = makeUniformMesh(N, M, type, theta)
bctype = [...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ClampedVel, ...
    enumBoundaryCondition.ClampedVel];

if (type == enumStdCell.Tri)
    mesh = makeUniformRotationTriMesh(N, [0, 1000], [-100, 100], ...
        M, ceil(M/5), bctype, 500, 0, theta);
elseif(type == enumStdCell.Quad)
    mesh = makeUniformRotationQuadMesh(N, [0, 1000], 0.2*[-100, 100], ...
        M, ceil(0.2*ceil(M/5)), bctype, 500, 0, theta );
else
    msgID = 'DamBreakDryUniformMesh:inputCellTypeError';
    msgtext = 'The input cell type should be NdgCellType.Tri or NdgCellType.Quad.';
    ME = MException(msgID, msgtext);
    throw(ME);
end
end% func
