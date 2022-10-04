classdef FlowOverATrapezoidTrech < SWEBarotropic3d
    %FLOWOVERATRAPEZOIDTRECH 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties(Constant)
        hcrit = 0.01;
        ChLength = 2.8
        ChWidth = 0.15
        startTime = 0
        finalTime = 100
        GotmFile = fullfile([pwd,'/Application/SWE/SWE3d/Benchmark/@FlowOverATrapezoidTrech/gotmturb.nml'])
    end
    
    properties
        SurfaceBoundaryEdgeType = 'Dirichlet'
    end    
    
    methods
        function obj = FlowOverATrapezoidTrech(N, Nz, M, Mz)
            %FLOWOVERATRAPEZOIDTRECH 构造此类的实例
            %   此处显示详细说明
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [ 1 2 4 6];
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
%             obj.varFieldIndex = [1 2 11];
            obj.varFieldIndex = [1 2];
%             obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.mesh3d );
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
                TempBottomElevation = zeros(mesh2d.cell.Np, mesh2d.K);
                Index = all(obj.meshUnion(1).mesh2d.x<=0.5);
                TempBottomElevation(:,Index) =  0.2;
                Index = all(obj.meshUnion(1).mesh2d.x>=0.5 & obj.meshUnion(1).mesh2d.x<=0.9);
                TempBottomElevation(:,Index) =  0.2 - (obj.meshUnion(1).mesh2d.x(:,Index) - 0.5) * 0.5;
                Index = all(obj.meshUnion(1).mesh2d.x>=0.9 & obj.meshUnion(1).mesh2d.x<=1.9);
                TempBottomElevation(:,Index) = 0;
                Index = all(obj.meshUnion(1).mesh2d.x>=1.9 & obj.meshUnion(1).mesh2d.x<=2.3);
                TempBottomElevation(:,Index) = (obj.meshUnion(1).mesh2d.x(:,Index) - 1.9) * 0.5;
                Index = all(obj.meshUnion(1).mesh2d.x>=2.3);
                TempBottomElevation(:,Index) = 0.2;
                fphys2d{m}(:, :, 4) = TempBottomElevation;
                %water depth
                fphys2d{m}(:,:,1) = 0.4 - fphys2d{m}(:, :, 4);
                %                  fphys2d{m}(:,:,1) = 2.89677;
            end
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
            % For 3d external field, the variable is organized as hu hv and
            % h.
            ks = 0.002;
            hu3d = zeros(size(obj.fext3d{1}(:,:,1)));
            [fm, ~] = obj.meshUnion.BoundaryEdge.matEvaluateSurfValue(fphys);
            Temph3d = fm(:,:,4);
%             Bz = Temph3d .* ( 1 + obj.meshUnion.BoundaryEdge.zb );
            Bz = 0.2 .* ( 1 + obj.meshUnion.BoundaryEdge.zb );
            Index = (Bz ~= 0 );
%             hu3d(Index) = Temph3d(Index) .* 0.033/0.4.*log(30*Bz(Index)/ks);
            hu3d(Index) = 0.2 .* 0.033/0.4.*log(30*Bz(Index)/ks);
            Index = any(Bz == 0);
            Tempz = Bz(:,Index)./2;
%             hu3d(1:2,Index) = 2*Temph3d(3:4,Index) .* 0.033/0.4.*log(30*Tempz(3:4,:)/ks) - ...
%                 hu3d(3:4,Index);
            hu3d(1:2,Index) = 2*0.2 .* 0.033/0.4.*log(30*Tempz(3:4,:)/ks) - ...
                hu3d(3:4,Index);
            obj.fext3d{1}(:,:,1) = hu3d;
            obj.fext2d{1}(:,:,1) = obj.meshUnion.BoundaryEdge.VerticalColumnIntegralField(hu3d);
            obj.fext3d{1}(:,:,3) = 0.2;
            obj.fext2d{1}(:,:,3) = 0.2;
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
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
            option('EosType') = enumEOSType.Jackett05;
            option('PhysicalSurfaceRoughnessLength') = 0.002;
            option('PhysicalBottomRoughnessLength') = 0.002;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.Constant;
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

end