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
        miu = 0.01
    end
    
    methods
        function obj = SecondOrderOperatorTest1d(N, M)
            [ obj.meshUnion ] = makeUniformMesh( N, M );
%             obj.meshUnion = obj.mesh3d;
            obj.Nmesh = 1;
            %             obj.outputFieldOrder2d = [];
            obj.fphys = obj.setInitialField;
            obj.option = obj.setOption( obj.option );
            obj.outputFile = obj.matInitOutput;
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
            fphys{1}(:,:,1) = 1/sqrt(4*0+1)*exp(-(obj.meshUnion(1).x-0.5).^2/obj.miu/(4*0+1));
%             fphys{1}(:,:,1) = 1/obj.miu*exp(-(obj.mesh3d(1).z+0.5).^2);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 2;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            option('EddyViscosityType') = enumEddyViscosity.Constant;
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

