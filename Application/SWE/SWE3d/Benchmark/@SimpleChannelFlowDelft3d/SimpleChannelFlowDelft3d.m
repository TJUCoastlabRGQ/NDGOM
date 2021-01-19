classdef SimpleChannelFlowDelft3d < SWEBarotropic3d
    %SIMPLECHANNELFLOWDELFT3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 10000;
        ChWidth = 2000;
        startTime = 0;
        %> final time

        finalTime = 7200;
        hcrit = 0.001;
        % to be corrected
        GotmFile = fullfile([pwd,'/Application/SWE/SWE3d/Benchmark/@SimpleChannelFlowDelft3d'],'/gotmturb.nml');        
    end
    
    methods
        function obj = SimpleChannelFlowDelft3d( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 3 4 5 10];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            obj.Cf{1} = 0.01 * ones(size(obj.Cf{1}));
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
                fphys2d{m}(:, :, 4) =  -0.0001 * mesh2d.x;
                %water depth
                fphys2d{m}(:,:,1) = 3.89677 + 0.025 - 0.0001 * (10000-mesh2d.x);
%                  fphys2d{m}(:,:,1) = 2.89677;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
            % For 3d external field, the variable is organized as hu hv and
            % h.
%             obj.BotBoundNewmannDate(:,:,1)  = obj. 
%            VCV = obj.meshUnion(1).cell.VCV;
%            Nz = obj.meshUnion(1).Nz;
%            Hu = VCV * fphys{1}(:,Nz:Nz:end,1);
%            Hv = VCV * fphys{1}(:,Nz:Nz:end,2);
%            H  = VCV * fphys{1}(:,Nz:Nz:end,4);
%            obj.BotBoundNewmannDate(:,:,1) = obj.Cf{1} .* sqrt( (Hu./H).^2 + ...
%                (Hv./H).^2 ) .* ( Hu./H ) * (-1);
%            obj.BotBoundNewmannDate(:,:,2) = obj.Cf{1} .* sqrt( (Hu./H).^2 + ...
%                (Hv./H).^2 ) .* ( Hv./H ) * (-1);
           
           hu3d = zeros(size(obj.fext3d{1}(:,:,1)));
           hu2d = zeros(size(obj.fext2d{1}(:,:,1)));
           h3d = zeros(size(obj.fext3d{1}(:,:,1)));
           h2d = zeros(size(obj.fext2d{1}(:,:,1)));
           Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedVel );
           hu3d(:,Index) = 5;
           obj.fext3d{1}(:,:,1) = hu3d;
           Index = ( obj.meshUnion(1).BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
           h3d(:,Index) = 3.89677 + 0.025;
           obj.fext3d{1}(:,:,3) = h3d;
           
           Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedVel );
           hu2d(:,Index) = 5;
           obj.fext2d{1}(:,:,1) = hu2d;
           Index = ( obj.mesh2d.BoundaryEdge.ftype == enumBoundaryCondition.ClampedDepth );
           h2d(:,Index) = 3.89677 + 0.025;
           obj.fext2d{1}(:,:,3) = h2d;
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 86400;
            outputIntervalNum = 3000;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantVerticalEddyViscosityValue')  = 0;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 0;
        end
        
    end    
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )


bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.ClampedVel, ...
    enumBoundaryCondition.ClampedDepth ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength], [ 0, obj.ChWidth ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end

