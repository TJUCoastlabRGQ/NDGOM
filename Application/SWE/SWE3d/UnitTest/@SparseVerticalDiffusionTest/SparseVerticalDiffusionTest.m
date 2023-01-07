classdef SparseVerticalDiffusionTest < SWEBarotropic3d
    %FLOWOVERATRAPEZOIDTRECH 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        hcrit = 0.01;
        ChLength = 10
        ChWidth = 10
        startTime = 0
        finalTime = 1
    end
    
    properties
        miu = 0.05
        
        Cexact
        
        DiffCexact
        
        DirichExact
        
        NewmannExact
    end
    
    
    properties
        SurfaceBoundaryEdgeType = 'Dirichlet'
    end
    
    methods
        function obj = SparseVerticalDiffusionTest(N, Nz, M, Mz)
            %FLOWOVERATRAPEZOIDTRECH 构造此类的实例
            %   此处显示详细说明
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 11];
            obj.Prantl = [1.0, 1.0, 1.0];
            obj.Nvar = 3;
            obj.varFieldIndex = [1, 2, 11];

            obj.matGetFunction;
            
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.mesh3d );
            
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            obj.fphys2d{1}(:,:,5) = ones(size(obj.fphys2d{1}(:,:,5)));
        end
        
        matTimeSteppingLai( obj );
        
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
                t = 0;
                z = mesh3d.z;
                fphys{m}(:,:,1) = eval(obj.Cexact);
                fphys{m}(:,:,11) = eval(obj.Cexact);
                fphys2d{m}(:, :, 4) = -1;
                %water depth
                fphys2d{m}(:,:,1) = 0 - fphys2d{m}(:, :, 4);
            end
        end
        
        function matUpdateExternalField( obj, time, ~, ~ )
            %Top first column, bottom second column
            t = time;
            z = zeros(obj.mesh2d.cell.Np,obj.mesh2d.K);
            obj.SurfBoundNewmannDate(:,:,1) = obj.miu * eval(obj.DiffCexact);
            obj.SurfBoundNewmannDate(:,:,3) = obj.miu * eval(obj.DiffCexact);
            %                         obj.NewmannExact(1) = obj.miu*1;
            z = -1 * ones(obj.mesh2d.cell.Np,obj.mesh2d.K);
            obj.BotBoundNewmannDate(:,:,1) = -1.0 * obj.miu * eval(obj.DiffCexact);
            
%             obj.BotBoundNewmannDate(:,:,3) = -1.0 * obj.miu * eval(obj.DiffCexact);
        end
        
        function matGetFunction(obj)
            syms z t;
            z0 = -0.5;
            obj.Cexact = 1/sqrt(4*t+1)*exp(-(z-z0)^2/obj.miu/(4*t+1));
            obj.DiffCexact = diff(obj.Cexact, z);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = obj.finalTime;
            outputIntervalNum = 1500;
            option('startTime') = obj.startTime;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
%             option('GOTMSetupFile') = obj.GotmFile;
            option('EosType') = enumEOSType.Jackett05;
            option('PhysicalSurfaceRoughnessLength') = 0.002;
            option('PhysicalBottomRoughnessLength') = 0.002;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = obj.miu;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 0.1;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ClampedVel, ...
    enumBoundaryCondition.ClampedDepth ];

% bctype = [ ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh(N, [0, obj.ChLength],...
    [-obj.ChWidth/2, obj.ChWidth/2], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

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

