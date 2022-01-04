classdef EkmanSpiral < SWEBarotropic3d
    %EKMANSPIRAL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 40000000;
        
        ChWidth  = 40000000;
        %> channel depth
        H0 = 200;
        %> x range
        %> start time
        startTime = 0;
        %> final time
        finalTime = 86400 * 100;
        
        hcrit = 0.001;
    end
    
    methods
        function obj = EkmanSpiral(N, Nz, M, Mz)
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 ];
            obj.outputFieldOrder3d = [ 1 2 4];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            obj.Cf{1} = 0.0025*ones(size(obj.mesh2d(1).x));
            
            obj.SurfBoundNewmannDate(:,:,2) = 0.1/1000 * ones(size(obj.SurfBoundNewmannDate(:,:,1)));%0.1
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
%                 Index = ( mesh3d.z == 0 );
%                 fphys{m}(Index) = 1.5;
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0 .* ones(mesh2d.cell.Np, mesh2d.K);
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
%            VCV = obj.meshUnion(1).cell.VCV;
%            Nz = obj.meshUnion(1).Nz;
%            Hu = VCV * fphys{1}(:,Nz:Nz:end,1);
%            Hv = VCV * fphys{1}(:,Nz:Nz:end,2);
%            H  = VCV * fphys{1}(:,Nz:Nz:end,4);
%            obj.BotBoundNewmannDate(:,:,1) = obj.Cf{1} .* sqrt( (Hu./H).^2 + ...
%                (Hv./H).^2 ) .* ( Hu./H ) * (-1);
%            obj.BotBoundNewmannDate(:,:,2) = obj.Cf{1} .* sqrt( (Hu./H).^2 + ...
%                (Hv./H).^2 ) .* ( Hv./H ) * (-1);           
        end        
        
        function [ option ] = setOption( obj, option )
            ftime = obj.finalTime;
            outputIntervalNum = 10000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.1;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 0;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.Beta;
            option('f0 for beta coriolis solver') = 2*sin(pi/4)*7.29*10^(-5);
            option('beta for beta coriolis solver') = 0.0;
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
end

