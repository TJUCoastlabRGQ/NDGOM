classdef EllipticProblem3d < SWEBarotropic3d
    %STANDINGWAVEINACLOSECHANNEL �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties ( Constant )
        %> channel length
        hcrit = 0.01;
%         ChLength = 100;
        ChLength = 2;
        %> channel width
        ChWidth = 2;
        %> channel depth
        H0 = 10;
        %> x range
        %> start time
        startTime = 0;
%         %> final time
        finalTime = 10;
    end
    
    properties
        dt
        miu0
        Lambda = 20;
        A = 0.1;
    end
    
    methods
        function obj = EllipticProblem3d( N, Nz, M, Mz )
            % setup mesh domain
            [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [1 2 3 11];
            obj.Nfield = 11;
            obj.Nvar = 3;
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, mesh3d );
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( mesh2d, mesh3d );
            %> time interval
            obj.dt = 0.005;
%             obj.Cf{1} = 0.0025/1000;
            obj.Cf{1} = 0*ones(size(mesh2d.x));
        end
        
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;                
                %water depth
                fphys2d{m}(:,:,1) =  obj.A * cos(2*pi*obj.mesh2d(m).x/obj.Lambda) - fphys2d{m}(:, :, 4);
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
            option('outputNcfileNum') = 20;                  
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
%             option('EddyViscosityType') = enumEddyViscosity.Constant;
%             option('GOTMSetupFile') = obj.GotmFile;
%             option('equationType') = enumDiscreteEquation.Strong;
%             option('integralType') = enumDiscreteIntegral.GaussQuadrature;
%             option('outputType') = enumOutputFile.VTK;
%             option('ConstantEddyViscosityValue') = 0;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 1; 
            
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
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
end
