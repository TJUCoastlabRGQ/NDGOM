classdef BottomBoundaryLayerCase < SWEBarotropic3d
    %STANDINGWAVEINACLOSECHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 10000;
        %> channel depth
        H0 = 15;
        %> x range
        %> start time
        startTime = 0;
        %> final time
        finalTime = 86400;
        hcrit = 0.001;
        % to be corrected
        GotmFile = fullfile([pwd,'/Application/SWE/SWE3d/Benchmark/@BottomBoundaryLayerCase/gotmturb.nml']);
    end
    
    properties
        dt
    end
    
    methods
        function obj = BottomBoundaryLayerCase( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3];
            obj.outputFieldOrder3d = [ 1 2 3 4 5];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
            obj.dt = 0.6;
            obj.Cf{1} = 0*ones(size(obj.mesh2d.x));
            %             obj.Cf{1} = 0.0025/1000;
            %             obj.Cf{1} = 0.0025;
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
                fphys2d{m}(:,:,1) = -mesh2d.x *10^(-5) - fphys2d{m}(:, :, 4);
                %                 Index = ( mesh2d.x <= 0 );
                %                 fphys2d{m}(Index) = 0.15/5000 * (mesh2d.x(Index) + 5000)./5000;
                %                 fphys2d{m}(~Index) = -0.15/5000 * (mesh2d.x(~Index) - 5000)./5000;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
           
           h3d = zeros(size(obj.fext3d{1}(:,:,1)));
           h2d = zeros(size(obj.fext2d{1}(:,:,1)));
           Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
           h3d(:,Index) = -obj.meshUnion(1).BoundaryEdge.xb(:,Index) *10^(-5)  + 15;
           obj.fext3d{1}(:,:,3) = h3d;
           
           Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
           h2d(:,Index) = -obj.meshUnion(1).mesh2d.BoundaryEdge.xb(:,Index) *10^(-5)  + 15;
           obj.fext2d{1}(:,:,3) = h2d;
        end        
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 7500;
            option('startTime') = 0.0;
            option('finalTime') = obj.finalTime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = obj.finalTime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            %             option('EddyViscosityType') = enumEddyViscosity.GOTM;
            %             option('GOTMSetupFile') = obj.GotmFile;
            %             option('equationType') = enumDiscreteEquation.Strong;
            %             option('integralType') = enumDiscreteIntegral.QuadratureFree;
            %             option('outputType') = enumOutputFile.VTK;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.01;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 0;
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.ClampedDepth, ...
    enumBoundaryCondition.ClampedDepth, ...
    enumBoundaryCondition.ClampedDepth, ...
    enumBoundaryCondition.ClampedDepth ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -obj.ChLength/2, obj.ChLength/2 ], 0.02*[ -obj.ChLength/2, obj.ChLength/2 ], ceil(obj.ChLength/M), 0.02*ceil(obj.ChLength/M), bctype);

cell = StdPrismQuad( N, Nz );
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

