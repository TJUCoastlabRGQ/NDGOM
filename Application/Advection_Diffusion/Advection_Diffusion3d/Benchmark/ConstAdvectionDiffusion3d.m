classdef ConstAdvectionDiffusion3d < Adv_DiffAbstract3d
    %CONSTADVECTIONDIFFUSION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    properties
        M
        Mz
        N
        Nz
    end
    
    methods
        function obj = ConstAdvectionDiffusion3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.M = M;
            obj.Mz = Mz;
            obj.N = N;
            obj.Nz = Nz;
            obj.miu = 0.001;
            obj.u0 = 1;
            obj.v0 = 1;
            obj.w0 = 1;
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );            
        end
        
        function fext = matGetExtFunc(obj, time)
            fext = obj.getExtFunc(obj.meshUnion(1), time);
        end
    end

    methods ( Access = protected )
        
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh3d = obj.mesh3d(m);
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                fphys{m}(:,:,1) = obj.getExtFunc( mesh3d, 0);
                fphys{m}(:,:,2) = obj.u0 * ones(size(mesh3d.x));
                fphys{m}(:,:,3) = obj.v0 * ones(size(mesh3d.x));
                fphys{m}(:,:,4) = obj.w0 * ones(size(mesh3d.x));
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
%             obj.BoundaryEdgefp3d{1} = sin(2*pi*time)*sin(2*pi*obj.mesh3d.BoundaryEdge.xb).*...
%                 sin(pi*obj.mesh3d.BoundaryEdge.yb).*sin(2*pi*obj.mesh3d.BoundaryEdge.zb);
            obj.fext{1} = sin(2*pi*time)*sin(2*pi*obj.mesh3d.BoundaryEdge.xb).*...
                sin(pi*obj.mesh3d.BoundaryEdge.yb).*sin(2*pi*obj.mesh3d.BoundaryEdge.zb);
%             obj.SurfaceBoundaryEdgefp3d{1} = sin(2*pi*time)*sin(2*pi*obj.mesh2d.x).*...
%                 sin(pi*obj.mesh2d.y).*sin(2*pi*1 );
            obj.SurfEdgefext{1} = sin(2*pi*time)*sin(2*pi*obj.mesh2d.x).*...
                sin(pi*obj.mesh2d.y).*sin(2*pi*1 );
%             obj.BottomBoundaryEdgefp3d{1} = sin(2*pi*time)*sin(2*pi*obj.mesh2d.x).*...
%                 sin(pi*obj.mesh2d.y).*sin(2*pi*0);
            obj.BotEdgefext{1} = sin(2*pi*time)*sin(2*pi*obj.mesh2d.x).*...
                sin(pi*obj.mesh2d.y).*sin(2*pi*0);
            obj.SurfBoundNewmannDate(:,:,1) = 2 * pi * obj.miu * sin(2*pi*time)*sin(2*pi*obj.mesh2d.x).*...
                sin(pi*obj.mesh2d.y).*cos(2*pi*1) .* 1;
            obj.BotBoundNewmannDate(:,:,1) = 2 * pi * obj.miu * sin(2*pi*time)*sin(2*pi*obj.mesh2d.x).*...
                sin(pi*obj.mesh2d.y).*cos(2*pi*0) .* (-1);            
%         BotBoundNewmannDate
        end
        
        function matEvaluateSourceTerm( obj, time )
            obj.frhs{1} = obj.frhs{1} + ...
                2*pi* (cos(2*pi*time).*sin(2*pi*obj.meshUnion.x) + cos(2*pi*obj.meshUnion.x).*sin(2*pi*time)).* sin(pi*obj.meshUnion.y).*sin(2*pi*obj.meshUnion.z)...
                + pi*(cos(pi*obj.meshUnion.y).*sin(2*pi*obj.meshUnion.z)+2*cos(2*pi*obj.meshUnion.z).*sin(pi*obj.meshUnion.y))*sin(2*pi*time).*sin(2*pi*obj.meshUnion.x) + ...
                9*obj.miu*pi^2*sin(2*pi*time).*sin(2*pi*obj.meshUnion.x).*sin(pi*obj.meshUnion.y).*sin(2*pi*obj.meshUnion.z);
        end
        
        function f_ext = getExtFunc( obj, mesh, time )
            f_ext = sin(2*pi*time)*sin(2*pi*mesh.x).*sin(pi*mesh.y).*...
                    sin(2*pi*mesh.z);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 0.2;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            option('AdvDiffVerticalDiffusionType') = enumVerticalDiffusion.Constant;
            dx = (obj.mesh2d.cell.r(2) - obj.mesh2d.cell.r(1))/2*(2/obj.M);
            dthu = min( 1/(2*obj.N+1) *  dx/obj.u0, 1/(2*obj.Nz+1) * dx^2/obj.miu);
            dthv = min( 1/(2*obj.N+1) *  dx/obj.v0, 1/(2*obj.Nz+1) * dx^2/obj.miu); 
            dtz = min( 1/(2*obj.Nz+1) *  dx/obj.w0, 1/(2*obj.Nz+1) * dx^2/obj.miu); 
            option('timeInterval') = min(min(dthu, dthv), dtz);
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
zs = ones(mesh2d.Nv, 1); zb = -1 * ones(mesh2d.Nv, 1);
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end