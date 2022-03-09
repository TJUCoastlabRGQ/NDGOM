classdef CouetteFlow < SWEBarotropic3d
    %STANDINGWAVEINACLOSECHANNEL �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties ( Constant )
        %> channel length
        ChLength = 2000;
        %> channel depth
        H0 = 20;
        %> x range
        %> start time
        startTime = 0;
        %> final time
        finalTime = 86400;
        hcrit = 0.001;
        % to be corrected
        GotmFile = fullfile([pwd,'/Application/SWE/SWE3d/Benchmark/@CouetteFlow/gotmturb.nml']);
    end
    
    methods
        function obj = CouetteFlow( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.Nfield = 15; % [ Hu, Hv, w, H, miu, Z, eta, Zx, Zy ]
            obj.fieldName3d = {'hu','hv','omega', 'h','nv','z','eta','zx','zy','w', 'hw','hc','rho', 'Tke', 'Eps'};
            obj.outputFieldOrder2d = [ 1 2 3];
            obj.outputFieldOrder3d = [ 1 2 3 4 5 14 15 ];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
            obj.Cf{1} = 0*ones(size(obj.mesh2d.x));
            
            obj.SurfBoundNewmannDate(:,:,1) = 1.027/1027 * ones(size(obj.mesh2d.x));%0.1
            %             obj.Cf{1} = 0.0025/1000;
            %             obj.Cf{1} = 0.0025;
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
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
           
        end        
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 750;
            option('startTime') = 0.0;
            option('finalTime') = obj.finalTime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = obj.finalTime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
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
            option('PhysicalSurfaceRoughnessLength') = 0.0030000000260770321;
            option('PhysicalBottomRoughnessLength') = 0.0030000000260770321;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
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
    [ -obj.ChLength/2, obj.ChLength/2 ], 1*[ -obj.ChLength/2, obj.ChLength/2 ], ceil(obj.ChLength/M), 1*ceil(obj.ChLength/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
[ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
[ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end

