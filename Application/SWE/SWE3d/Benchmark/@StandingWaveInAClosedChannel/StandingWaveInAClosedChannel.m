classdef StandingWaveInAClosedChannel < SWEBarotropic3d
    %STANDINGWAVEINACLOSECHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        hmin = 0.001;
        %> channel length
        ChLength = 100;
        %> channel width
        ChWidth = 20;
        %> channel depth
        H0 = 10;
        %> x range
        %> start time
        startTime = 0;
        %> final time
        finalTime = 30;
        %> casename
        casename = 'StandingWaveInAClosedChannel';
        Cf = 1e2;
    end
    
    properties
        dt
        Taux
        Tauy
    end
    
    methods
        function obj = StandingWaveInAClosedChannel( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            % set initilize physical field
            [ obj.fphys2d, obj.fphys3d ] = SetInitialField( obj );
            %> time interval
            obj.dt = 0.001;
            
            obj.miu = 10^(-2);
            
            obj.Taux = zeros(size(obj.fphys2d{1}(:,:,1)));
            obj.Tauy = zeros(size(obj.fphys2d{1}(:,:,1)));
        end
        
        AnalysisResult2d( obj );
        AnalysisResult3d( obj );
        
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [fphys2d, fphys3d] = SetInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys3d = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys3d{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield3d );
                
                Lambda = 20;
                % surface elevation
                fphys2d{m}(:,:,1) =  0.1 * cos(2*pi*mesh2d.x/Lambda);
                % bottom elevation
                fphys2d{m}(:, :, 5) = -obj.H0;
                % water depth
                fphys2d{m}(:, :, 4) = fphys2d{m}(:,:,1) - fphys2d{m}(:, :, 5);
                % water depth
                fphys3d{m}(:, :, 6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
            end
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
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, obj.ChWidth/(obj.ChLength/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end

