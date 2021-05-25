classdef PassiveTransportTest3d < SWEBarotropic3d
    %PASSIVETRANSPORTTEST3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        Depth = 1.0
        u = 0.5
        FinalTime = 10
    end
    
    properties
        Cexact
    end
    
    properties
        ChLength = 30
        ChWidth = 2
    end
    
    properties(Constant)
        hcrit = 0.05
    end
    
    methods
        function obj = PassiveTransportTest3d(N, Nz, M, Mz)
            [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.outputFieldOrder3d = [1 2 12];
            obj.Nvar = 3;
            obj.varFieldIndex = [1 2 12];            
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.matGetExactData;
        end
    end
    
    methods(Access = protected)
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            syms x;
            C = sech(x);
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                fphys{m}(:,:,1) = obj.Depth*obj.u.*ones(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K);
                x = obj.meshUnion(m).x;
                fphys{m}(:,:,12) = obj.Depth*eval(C);
                fphys2d{m}(:,:,1) = obj.Depth*ones( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K);
            end
        end
        
        function matGetExactData(obj)
            syms x t;
            C = sech(x);
            obj.Cexact = zeros(size(obj.meshUnion.x));
            mesh = obj.meshUnion;
            xmax = max( max(mesh.x) );
            xmin = min( min(mesh.x) );
            
            T = mod(obj.FinalTime, (xmax - xmin)/obj.u);
            Index = mesh.x - xmin >= obj.u * T;
            obj.Cexact(Index) = obj.Depth * double(subs(C,{x,t},{ mesh.x(Index) - obj.u * T,  0}));
            
            Index = mesh.x - xmin <= obj.u * T;
            obj.Cexact(Index) = obj.Depth * double(subs(C,{x,t},{ mesh.x(Index) - obj.u * T + xmax - xmin,  0}));
        end
        
        function [ option ] = setOption( obj, option )
            ftime = obj.FinalTime;
            outputIntervalNum = 2500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.01;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 0.1;
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
    [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

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