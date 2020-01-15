classdef TidalWave< SWEBarotropic3d
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
        hmin = 0.001;
    end
    
    properties
        dt
    end
    
    methods
        function obj = TidalWave( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3];
            obj.outputFieldOrder = [ 1 2 3 4];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            obj.WindTaux{1} = zeros(size(obj.mesh2d(1).y));
            obj.WindTauy{1} = zeros(size(obj.mesh2d(1).y));
            obj.Cf{1} = 0*ones(size(obj.mesh2d(1).x));
        end
        
        EntropyAndEnergyCalculation(obj);
        
        AnalysisResult2d( obj );
        AnalysisResult3d( obj );
        
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0*ones(mesh2d.cell.Np, mesh2d.K);
            end
        end
        
        function matUpdateExternalField( obj, time, ~, ~ )
            A = 0.5;
            omega = 2*pi/12/3600;
            L = 7.4110e+05;
            k = 2*pi/L;
            for m = 1:obj.Nmesh
                obj.fext2d{m}(:, :, 1) = obj.H0 + A*cos(omega*time-k*L + pi/2);
                obj.fext3d{m}(:, :, 4) = obj.H0 + A*cos(omega*time-k*L + pi/2);
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
            option('temporalDiscreteType') = enumTemporalDiscrete.ARK232;
            option('EddyViscosityType') = enumEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantEddyViscosityValue') = 0;
        end
        
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

% bctype = [ ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.ZeroGrad, ...
%     enumBoundaryCondition.ZeroGrad ];

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ClampedDepth ];

mesh2d = makeUniformTriMesh( N, ...
    [ 0, obj.ChLength ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

