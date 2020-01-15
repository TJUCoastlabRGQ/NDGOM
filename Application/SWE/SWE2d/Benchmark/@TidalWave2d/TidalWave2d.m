classdef TidalWave2d< SWEConventional2d
    %WINDDRIVENFLOW 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 24000;
        ChWidth = 1200;
        %> channel depth
        H0 = 30;
        %> start time
        startTime = 0;
        %> final time
        finalTime = 216000;
    end
    
    properties
        dt
    end
    
    methods
        function obj = TidalWave2d( N, M )
            [ mesh ] = makeUniformMesh(N, M);
            obj = obj@SWEConventional2d();
            obj.hmin = 1e-3;      
            obj.initPhysFromOptions( mesh );
        end
        
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function fphys = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = zeros( mesh.cell.Np, mesh.K, obj.Nfield );
                % bottom elevation
                fphys{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys{m}(:, :, 1) = obj.H0*ones(mesh.cell.Np, mesh.K);
            end
        end
        
        function matUpdateExternalField( obj, time, ~, ~ )
            A = 0.5;
            omega = 2*pi/12/3600;
            L = 7.4110e+05;
            k = 2*pi/L;
            for m = 1:obj.Nmesh
                obj.fext{m}(:, :, 1) = obj.H0 + A*cos(omega*time-k*L + pi/2);
            end
            
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 216000;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
        end
        
    end
    
end

function mesh = makeUniformMesh( N, M)

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ClampedDepth ];

mesh = makeUniformTriMesh( N, [ 0, 24000 ], [ -600, 600 ], ceil(24000/M), ceil(1200/M), bctype);
end

