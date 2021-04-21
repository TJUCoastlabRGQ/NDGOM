classdef Limiter3dTest < SWEBarotropic3d
    %LIMITER3DTEST 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 0.5;
        ChWidth = 0.5;
        %> channel depth
        H0 = 1;
        %> start time
        startTime = 0;
        %> final time

        finalTime = 3600;
        hcrit = 0.001;
    end
    
    properties
        dt
        
        VertSolver
        
        Limiter
    end
    
    methods
        function obj = Limiter3dTest( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.VertSolver = SWE3dVerticalVelocitySolver( obj.mesh2d, obj.mesh3d );
            obj.Limiter =  NdgVertLimiter3d( obj.mesh3d );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 3 10];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
        end
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
                fphys{m} = 0.1 * ones( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                fphys{m}(1,1,1) = 0.05;
                fphys{m}(4,1,1) = 0.15;
                fphys{m}(1,4,1) = 0.05;
                fphys{m}(4,4,1) = 0.15;
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0 .* ones(mesh2d.cell.Np, mesh2d.K);
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
         
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 3600;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantVerticalEddyViscosityValue') = 0.01;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 100;
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
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformTriMesh( N, ...
    [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end
