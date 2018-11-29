classdef OpenChannel3d < LSWEAbstract3d
    %OPENCHANNEL3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        %> channel length
        ChLength = 100e3;
        %> channel width
        ChWidth = 20e3;
        %> channel depth
%         H = 12;
        H = 5;
        % frequency
        omega = 2 * pi / 12.4 / 3600;
        %> start time
        startTime = 0;
        %> final time
        finalTime = 24 * 3600;
        %> output time interval
        outputTimeInterval = 240;
        %> casename
        casename = 'TidalForcingWithHorizontalBottomRectangularBasin';
        %> amplitude
        a = 0.25;
%         miu0 = 1e-4; 
        miu0 = 0.0015/5/5; 
    end
    
    methods
        function obj = OpenChannel3d( N, Nz, M, Mz )
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            % initilize physical field
            [ obj.fphys2d, obj.fphys3d ] = obj.setInitialField( );
            [ obj.fext2d, obj.fext3d ] = obj.initBoundaryField( );
            % set vertical viscosity
%             obj.miu = 0;
            obj.miu = 0.0015;
            %> linear slip parameter
%             obj.K = 0;
            obj.K = 0.01;
            %> time interval
            obj.dt = 0.000001;
        end
        
        AnalysisResult2d( obj );
        AnalysisResult3d( obj );
    end
    
    methods ( Access = protected )
        %> evaluate boundary flux for 2d PCE
        frhs2d = EvaluatePCE2d_BoundaryKernel( obj, edge, fphys2d, fext2d )
        %> evaluate horizontal flux for 3d momentum equations
        frhs3d = Evaluate3d_HorizontalBoundaryKernel( obj, edge, fphys3d, fext3d )
        %> update external fields
        matUpdateExternalField( obj, time, fphys2d, fphys3d );
        
        function [ fext2d, fext3d ] = initBoundaryField( obj )
            fext2d = cell( obj.Nmesh, 1 );
            fext3d = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh2d = obj.mesh2d(m);
                edge2d = mesh2d.BoundaryEdge;
                fext2d{m} = zeros( edge2d.Nfp, edge2d.Ne, obj.Nfield2d );
                
                edge3d = obj.mesh3d(m).BoundaryEdge;
                fext3d{m} = zeros( edge3d.Nfp, edge3d.Ne, obj.Nfield3d );
                % fext3d{m}(:, :, 1) = fext3d{m}(:, :, 5) .^2 ./ obj.miu / 2 ...
                %    .* obj.phi / obj.ChLength .* ( edge3d.zb + 1 ) .^2;
            end
        end
        
        function [ fphys2d, fphys3d ] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys3d = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.mesh2d(m);
                mesh3d = obj.mesh3d(m);
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys3d{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield3d );
                
                % surface elevation
                fphys2d{m}(:, :, 1) = obj.a * ( ...
                    cos( obj.omega / sqrt( obj.gra * obj.H ) * ( obj.ChLength - mesh2d.x ) ) / ...
                    cos( obj.omega / sqrt( obj.gra * obj.H ) * obj.ChLength ) );
                % bottom elevation
                fphys2d{m}(:, :, 5) = - obj.H;
                % water depth
                fphys2d{m}(:, :, 4) = obj.H;
                
                % water depth
                fphys3d{m}(:, :, 6) = mesh3d.Extend2dField( fphys2d{m}(:, :, 4) );
                % bottom elevation
                fphys3d{m}(:, :, 7) = mesh3d.Extend2dField( fphys2d{m}(:, :, 1) );
                
            end
        end
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformTriMesh( N, ...
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, 2, bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end
