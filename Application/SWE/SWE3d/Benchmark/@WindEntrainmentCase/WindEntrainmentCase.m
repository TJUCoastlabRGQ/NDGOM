classdef WindEntrainmentCase < SWEBaroclinic3d
    %WINDENTRAINMENTCASE 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        ChLength = 1000
        
        ChWidth = 1000
        
        H0 = 50
        
        finalTime = 1*3600
        
        GotmFile = fullfile([pwd,'/Application/SWE/SWE3d/Benchmark/@WindEntrainmentCase/gotmturb.nml'])
    end
    
    properties
        MLDfid
    end
    
    properties ( Constant )
        hcrit = 0.001;
    end
    
    methods
        function obj = WindEntrainmentCase(N, Nz, M, Mz)
            [ obj.mesh2d, obj.mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            
            obj.initPhysFromOptions( obj.mesh2d, obj.mesh3d );
            
            obj.SurfBoundNewmannDate(:,:,1) = 0.1027/obj.rho0 * ones(size(obj.mesh2d.x));%0.1
            
            obj.MLDfid = fopen('Result\WindEntrainmentCase\MixedLayerDepth.dat','w');
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
                %> For salinity
                fphys{m}(:,:,15) = obj.H0*obj.S0;
                %> For temperature
%                 Ttop = 20;
%                 fphys{m}(:,:,14) = obj.H0 * (Ttop - (0.01)^2*obj.H0/obj.alphaT/obj.gra/obj.rho0*(0-mesh3d.z));
                fphys{m}(:,:,14) = obj.H0 * obj.T0;
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;
                %water depth
                fphys2d{m}(:,:,1) = obj.H0;
            end
        end
        
        function matUpdateOutputResult( obj, time, fphys2d, fphys3d )
            matUpdateOutputResult@SWEAbstract3d( obj, time, fphys2d, fphys3d );
            Index = fphys3d{1}(:,:,14) > 1e-5;
            depth = -1 * obj.H0 * min(min(obj.meshUnion.z(Index)));
            fprintf(obj.MLDfid, '%12.8f  %12.8f\n', time, depth );
        end
        
        function matUpdateFinalResult( obj, time, fphys2d, fphys )
            matUpdateFinalResult@SWEAbstract3d( obj, time, fphys2d, fphys );
            fclose(obj.MLDfid);
        end
        
        function matUpdateExternalField( obj, time, fphys2d, fphys )
            %doing nothing
        end
        
        function [ option ] = setOption( obj, option )
            outputIntervalNum = 7500;
            option('startTime') = 0.0;
            option('finalTime') = obj.finalTime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = obj.finalTime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.GOTM;
            option('GOTMSetupFile') = obj.GotmFile;
            option('EosType') = enumEOSType.Linear;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.0001;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('PhysicalSurfaceRoughnessLength') = 0.003;
            option('PhysicalBottomRoughnessLength') = 0.003;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.None;
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
    [ -obj.ChLength/2, obj.ChLength/2 ], 1*[ -obj.ChLength/2, obj.ChLength/2 ], ceil(obj.ChLength/M), 1*ceil(obj.ChLength/M), bctype);

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