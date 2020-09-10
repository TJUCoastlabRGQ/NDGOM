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
        miu = 0.001
        
        Cexact
        
        DiffCexact
        
        DirichExact
        
        NewmannExact
    end
    
    methods
        function obj = SecondOrderOperatorTest1d(N, M)
            [ obj.meshUnion ] = makeUniformMesh( N, M );
            %             obj.meshUnion = obj.mesh3d;
            obj.Nmesh = 1;
            %             obj.outputFieldOrder2d = [];
            obj.matGetFunction;
            obj.option = obj.setOption( obj.option );
            obj.fphys = obj.setInitialField;
            obj.DirichExact = zeros(1,2);
            obj.NewmannExact = zeros(1,2);
            %             obj.outputFile = obj.matInitOutput;
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
            x = obj.meshUnion.x;
            t = 0;
            fphys{1}(:,:,1) = eval(obj.Cexact).*ones(size(obj.meshUnion.x));
            %             fphys{1}(:,:,1) = 1/obj.miu*exp(-(obj.mesh3d(1).z+0.5).^2);
        end
        
        function matGetFunction(obj)
            syms x t;
            x0 = 0.5;
%             obj.Cexact = 1/sqrt(4*t+1)*exp(-(x-x0)^2/obj.miu/(4*t+1));
%             obj.DiffCexact = diff(obj.Cexact, x);
%                         obj.Cexact = 0*t+0.*x;
            obj.Cexact =-1* sin(2*pi*x + 0 * t);
            %             obj.DiffCexact = diff(obj.Cexact, x);
            
        end
        
        function matUpdateExternalField( obj, time )
            t = time;
            x = 0;
%             obj.DirichExact(1) = eval(obj.Cexact);
%             obj.NewmannExact(1) = obj.miu * eval(obj.DiffCexact);
                        obj.NewmannExact(1) = obj.miu*1;
            x = 1;
%             obj.DirichExact(2) = eval(obj.Cexact);
%             obj.NewmannExact(2) = obj.miu * eval(obj.DiffCexact);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 2000;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            option('EddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('ConstantEddyViscosityValue') = 0.01;
        end
    end
    
end

function [ mesh ] = makeUniformMesh( N, M )
xlim = [0, 1];
bcType = [enumBoundaryCondition.SlipWall, enumBoundaryCondition.SlipWall];
[ mesh ] = makeUniformMesh1d( N, xlim, M, bcType );
end

