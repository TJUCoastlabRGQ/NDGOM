classdef ModeSplitStandingWaveInAClosedChannel < SWEMSBarotropic3d
    %STANDINGWAVEINACLOSECHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        hmin = 0.001;
        %> channel length
        ChLength = 100;
        %> channel width
        ChWidth = 20;
        %> channel depth
        H0 = 7.5;
        %> x range
        %> start time
        startTime = 0;
        %> final time
        finalTime = 20;
    end
    
    properties
        dt
%         Taux
%         Tauy
    end
    
    methods
        function obj = ModeSplitStandingWaveInAClosedChannel( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            
%             obj.Solver2d = ModeSplitStandingWaveInAClosedChannel2d( obj.mesh2d );
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            % set initilize physical field
            
            %> time interval
            obj.dt = 0.3;
            obj.outputFieldOrder2d = 4;
%             obj.miu = 0.001;
            obj.miu = 0.01;
            obj.Cf{1} = 1000;
            
            
        end
        
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
                
                Lambda = 20;
                % bottom elevation
                fphys2d{m}(:, :, 5) = -obj.H0;                
                % water depth
                fphys2d{m}(:,:,1) =  0.01 * cos(2*pi*mesh2d.x/Lambda) - fphys2d{m}(:, :, 5);
                % surface elevation
                fphys2d{m}(:, :, 4) = 0.01 * cos(2*pi*mesh2d.x/Lambda);
                % water depth
                fphys{m}(:, :, 6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 1) );
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
            option('outputNcfileNum') = 500;                  
            option('temporalDiscreteType') = enumTemporalDiscrete.Euler;
            option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
%             option('nonhydrostaticType') = enumNonhydrostaticType.Nonhydrostatic;            
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformTriMesh( N, ...
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, obj.ChWidth/(obj.ChLength/M), bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

