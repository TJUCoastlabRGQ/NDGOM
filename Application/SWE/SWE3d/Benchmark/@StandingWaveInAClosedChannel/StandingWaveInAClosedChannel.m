classdef StandingWaveInAClosedChannel < SWEBarotropic3d
    %STANDINGWAVEINACLOSECHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 100;
        %> channel width
        ChWidth = 2;
        %> channel depth
        H0 = 7.612;
        %> x range
        %> start time
        startTime = 0;
        %> final time
        finalTime = 65;
        hcrit = 0.001;
        % to be corrected
        GotmFile = fullfile('D:\PhdResearch\Application\SWE\SWE3d\Benchmark\@StandingWaveInAClosedChannel','\gotmturb.nml');
    end
    
    properties
        dt
        miu0
        Lambda = 20;
        A = 0.1;
    end
    
    methods
        function obj = StandingWaveInAClosedChannel( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            % set initilize physical field
%             [ obj.fphys2d, obj.fphys ] = obj.setInitialField;
            %> time interval
            obj.dt = 0.01;
            
%             obj.miu0{1} = 0.001;
            obj.miu0{1} = 0;
%             obj.miu = 0;
            obj.Cf{1} = 0.0025/1000;
%             obj.Cf{1} = 0;
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
                fphys2d{m}(:,:,1) =  obj.A * cos(2*pi*mesh2d.x/obj.Lambda) - fphys2d{m}(:, :, 4);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 65;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;                  
            option('temporalDiscreteType') = enumTemporalDiscrete.RK45;
            option('EddyViscosityType') = enumEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
%             option('limiterType') = enumLimiter.Vert;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, ceil(obj.ChWidth/(obj.ChLength/M)), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

