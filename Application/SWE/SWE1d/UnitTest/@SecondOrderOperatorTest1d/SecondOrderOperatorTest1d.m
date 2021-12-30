classdef SecondOrderOperatorTest1d < SWEConventional1d
    %SECONDORDEROPERATORTEST 此处显示有关此类的摘要
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
    end
    
    properties
        miu = 0.1
        
        Cexact
        
        DiffCexact
        
        DirichExact
        
        NewmannExact
        %The final exact solution of C at the final step
        FCexact
    end
    
    methods
        
        function obj = SecondOrderOperatorTest1d(N, M)
            obj = obj@SWEConventional1d();
            [ mesh ] = makeUniformMesh( N, M );        
            obj.matGetFunction;
            obj.fieldName = {'C'};
            obj.outputFieldOrder1d = 1;
            obj.initPhysFromOptions( mesh );
            obj.DirichExact = zeros(1,2);
            obj.NewmannExact = zeros(1,2);
            obj.matGetFinalSolution;
            obj.dt = 0.05;
        end        
        
        function drawComparison( obj )
            hold on;
            t = obj.getOption('finalTime');
            x = obj.meshUnion.x;
            ed = eval(obj.Cexact).*ones(size(obj.meshUnion.x));
            plot(x(1:5:end),ed(1:5:end),'ro','LineWidth',1.5);
%             str =['$\mu=$', num2str(obj.miu)];
            str = '$\frac{1}{\sqrt{(4t+1)}}e^{\frac{-(x-0.5)^2}{4t+1}}$';
            title(str,'Interpreter','Latex','fontsize',15);
            plot(x, obj.fphys{1}(:,:,1),'k','LineWidth',1.5);
            legend({'Exact','Simulation'});
            legend('boxoff');
%          title('$u=x(1-x)$','interpreter','Latex');
            box on;
            set(gca, 'Linewidth',1.5, 'Fontsize',12);            
        end
        
        matTimeStepping343(obj);
        
        matTimeStepping222(obj);
        
        matTimeExlicitStepping222(obj);
        
        matTimeExlicitSteppingForwardEuler(obj);
    end
    
    methods( Hidden )
        function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys )
            %an empty function
        end
    end
    
    methods(Access = protected)
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            fphys{1} = zeros(obj.meshUnion.cell.Np, obj.meshUnion.K, 3);
            x = obj.meshUnion.x;
            t = 0;
            fphys{1}(:,:,1) = eval(obj.Cexact).*ones(size(obj.meshUnion.x));
            %             fphys{1}(:,:,1) = 1/obj.miu*exp(-(obj.mesh3d(1).z+0.5).^2);
        end
        
        function matGetFinalSolution(obj)
            x = obj.meshUnion.x;
            t = obj.getOption('finalTime');
            obj.FCexact = eval(obj.Cexact).*ones(size(obj.meshUnion.x));            
        end
        
        function matGetFunction(obj)
            syms x t;
            x0 = 0.5;
            obj.Cexact = 1/sqrt(4*t+1)*exp(-(x-x0)^2/obj.miu/(4*t+1));
            obj.DiffCexact = diff(obj.Cexact, x);   
        end
        
        function matUpdateExternalField( obj, time )
            t = time;
            x = obj.meshUnion.x(1);
            obj.DirichExact(1) = eval(obj.Cexact);
            obj.NewmannExact(1) = obj.miu * eval(obj.DiffCexact);
%                         obj.NewmannExact(1) = obj.miu*1;
            x = obj.meshUnion.x(end);
            obj.DirichExact(2) = eval(obj.Cexact);
            obj.NewmannExact(2) = obj.miu * eval(obj.DiffCexact);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 50;
            outputIntervalNum = 5000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('EddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('ConstantEddyViscosityValue') = 0.01;
            option('limiterType') = enumLimiter.None;
            option('outputType') = enumOutputFile.NetCDF;
        end
    end
    
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 1];
bcType = [enumBoundaryCondition.SlipWall, enumBoundaryCondition.SlipWall];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end

