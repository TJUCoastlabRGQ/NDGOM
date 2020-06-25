classdef ConstAdvection3d < Adv_DiffAbstract3d
    %CONSTADVECTION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        x0
        y0
    end
    
    methods
        function obj = ConstAdvection3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.miu = 0;
            obj.u0 = 0.5;
            obj.v0 = 0.5;
            obj.x0 = -0.5;
            obj.y0 = -0.5;
            obj.w0 = 0;
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );            
        end
        
        function fext = matGetExtFunc(obj, time)
            fext = obj.getExtFunc(obj.meshUnion(1).x, obj.meshUnion(1).y, time);
        end
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m} = getExtFunc(obj, mesh.x, mesh.y, 0);
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            obj.BoundaryEdgefp3d{1} = obj.getExtFunc( obj.mesh3d.BoundaryEdge.xb, obj.mesh3d.BoundaryEdge.yb, time);
            obj.SurfaceBoundaryEdgefp3d{1} = obj.getExtFunc( obj.mesh2d.x, obj.mesh2d.y, time);
            obj.BottomBoundaryEdgefp3d{1} = obj.getExtFunc( obj.mesh2d.x, obj.mesh2d.y, time);
            obj.SurfBoundNewmannDate(:,:,1) =zeros(size(obj.mesh2d.x));  
            obj.BotBoundNewmannDate(:,:,1) = zeros(size(obj.mesh2d.x));            
%         BotBoundNewmannDate
        end
        
        function matEvaluateSourceTerm( obj, time )
                %doing nothing
        end
        
        function f_ext = getExtFunc( obj, x, y , time )
            xc = obj.x0 + obj.u0.*time;
            yc = obj.y0 + obj.v0.*time;
            
            sigma = 125*1e3/(33*33);
            t = -( (x-xc).^2+(y-yc).^2 )*sigma;
            f_ext = exp(t);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 2;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('AdvDiffVerticalDiffusionType') = enumVerticalDiffusion.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('AdvDiffConstantVerticalDiffusionValue') = obj.miu;
            option('AdvDiffHorizontalDiffusionType') = enumHorizontalDiffusion.Constant;
            option('AdvDiffConstantHorizontalDiffusionValue') = obj.miu;
        end
        
    end      
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -1, 1 ], [ -1, 1 ], M, M, bctype);

cell = StdPrismQuad( N, Nz );
zs = ones(mesh2d.Nv, 1); zb = -1 * zs;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end
