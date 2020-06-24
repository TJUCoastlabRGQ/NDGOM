classdef WindDrivenFlow < SWEBarotropic3d
    %WINDDRIVENFLOW 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        ChLength = 2000;
        ChWidth = 300;
        %> channel depth
        H0 = 10;
        %> start time
        startTime = 0;
        %> final time

        finalTime = 50;
        hcrit = 0.001;
    end
    
    properties
        dt
    end
    
    methods
        function obj = WindDrivenFlow( N, Nz, M, Mz )
            % setup mesh domain
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 3 10];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            %> time interval
            %             obj.dt = 0.02;
            obj.Cf{1} = 0.005*ones(size(obj.mesh2d(1).x));
            
            obj.SurfBoundNewmannDate(:,:,1) = 1.5/1000 * ones(size(obj.SurfBoundNewmannDate(:,:,1)));%0.1
            %             Index =( all( obj.mesh2d.x - obj.ChLength/2  + 2*M > -1e-5 ));
            %             obj.WindTaux{1}(:,Index) = 0;
            %             Index =( all(obj.mesh2d.x + obj.ChLength/2 - 2*M < 1e-5 ));
            %             obj.WindTaux{1}(:,Index) = 0;
            %             obj.Cf{1} = 0.0025/1000;
        end
        
        matExplicitRK222(obj);
        
        
        EntropyAndEnergyCalculation(obj);
        
        AnalysisResult2d( obj );
        AnalysisResult3d( obj );
        
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
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0 .* ones(mesh2d.cell.Np, mesh2d.K);
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
% %             obj.BotBoundNewmannDate(:,:,1)  = obj. 
           VCV = obj.meshUnion(1).cell.VCV;
           Nz = obj.meshUnion(1).Nz;
           Hu = VCV * fphys{1}(:,Nz:Nz:end,1);
           Hv = VCV * fphys{1}(:,Nz:Nz:end,2);
           H  = VCV * fphys{1}(:,Nz:Nz:end,4);
           obj.BotBoundNewmannDate(:,:,1) = obj.Cf{1} .* sqrt( (Hu./H).^2 + ...
               (Hv./H).^2 ) .* ( Hu./H ) * (-1);
           obj.BotBoundNewmannDate(:,:,2) = obj.Cf{1} .* sqrt( (Hu./H).^2 + ...
               (Hv./H).^2 ) .* ( Hv./H ) * (-1);           
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 50;
            outputIntervalNum = 100;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.Constant;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.Smagorinsky;
            option('ConstantHorizontalEddyViscosityValue') = 100;
        end
        
    end
    
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

% bctype = [ ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.SlipWall, ...
%     enumBoundaryCondition.ZeroGrad, ...
%     enumBoundaryCondition.ZeroGrad ];

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
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end

