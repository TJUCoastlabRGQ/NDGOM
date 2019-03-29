classdef LSWEStandingWaveInAClosedChannel < LSWEAbstract3d
      properties ( Constant )
        %> channel length
        ChLength = 100;
        %> channel width
        ChWidth = 20;
        %> channel depth
        H0 = 10;
        %> start time
        %> x range
        %> start time
        startTime = 0;
        %> final time
        finalTime = 30;
%         externTimeInterval = 5;
        %> output time interval
        outputTimeInterval = 240;
        %> casename
        casename = 'LSWEStandingWaveInAClosedChannel';
        %> amplitude
%         a = 0.25;
        cD = 1e2;
        miu0 = 1e-4;
    end
    
    properties
        timeExt
        etaExt
    end
    
    methods
        function obj = LSWEStandingWaveInAClosedChannel( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            % set initilize physical field
            [ obj.fphys2d, obj.fphys3d ] = SetInitialField( obj );
            % set parameter
            obj = SetParameter( obj );
            %> time interval
            obj.dt = 0.005;
            % setup external field
%             makeBoundaryField( obj );            
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
        [eta, u] = matEvaluateExactSolution( obj, xb, zb, time );
%         [ W ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d );
        [ dHUdx, dHVdy ] = obj.matEvaluatePartialDerivativeTermX( mesh3d,data );
%         function makeBoundaryField( obj )
%             % analytical boundary solution
%             Nstep = ceil( obj.finalTime / obj.externTimeInterval + 1 );
%             % edge2d = obj.mesh2d.BoundaryEdge;
%             % edge3d = obj.mesh3d.BoundaryEdge;
%             obj.etaExt = zeros( Nstep, 1 );
%             %uExt = zeros( edge3d.Nfp, edge3d.Ne );
%             obj.timeExt = linspace( 0, obj.finalTime, Nstep);
%             [ obj.etaExt, ~ ] = ...
%                 obj.matEvaluateExactSolution( obj.x2, 0, obj.timeExt );
%             % [ ~, obj.uExt ] = ...
%             %    obj.matEvaluateExactSolution( edge3d.xb, edge3d.zb, time(n) );
%             % obj.timeExt(n+1) = obj.timeExt(n) + obj.dt;
%         end
        
        function obj = SetParameter( obj )
            % set vertical viscosity
            obj.miu = obj.miu0 .* obj.fphys3d{1}(:, :, 6) .^2;
            % obj.miu = 0;
            %> linear slip parameter
            botEdge = obj.mesh3d.BottomBoundaryEdge;
%             ind = botEdge.FToN1 + (botEdge.FToE(1, :) - 1) .* obj.mesh3d.cell.Np;
            ind = bsxfun( @plus, botEdge.FToN1, (botEdge.FToE(1, :) - 1) .* obj.mesh3d.cell.Np);
            miu1 = obj.miu( ind );
            [ fm, ~ ] = botEdge.matEvaluateSurfValue( obj.fphys3d );
            
            obj.K = 0;
%             obj.K = obj.cD .*  miu1./ fm(:, :, 6);
        end
        
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

mesh2d = makeUniformTriMesh( N, ...
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, obj.ChWidth/(obj.ChLength/M), bctype);

cell = StdPrismTri( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1 );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1 );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end