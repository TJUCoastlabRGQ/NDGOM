classdef HorizontalSecondOrderTermTest < SWEBarotropic3d
    %STANDINGWAVEINACLOSECHANNEL 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties ( Constant )
        %> channel length
        hcrit = 0.01;
%         ChLength = 100;
        ChLength = 20;
        %> channel width
        ChWidth = 20;
        %> channel depth
        H0 = 10;
        %> x range
        %> start time
        startTime = 0;
%         %> final time
        finalTime = 10;
    end
    
    
    properties
        dt
        miu0
        Lambda = 20;
        A = 0.1;
    end
    
    properties
        
        UE
        
        RHSE
        
        ExactResult
        
        SimulationResult
    end
    
    methods
        function obj = HorizontalSecondOrderTermTest

           [ mesh2d, mesh3d ] = makeChannelMesh( obj, 1, 1, 100, 1 );
           
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, mesh3d );
            
            obj.initPhysFromOptions( mesh2d, mesh3d );
            
            obj.matGetFunction;

        end
        
        function EquationSolve(obj)
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            obj.ExactResult = eval(obj.UE);
            rhs = eval(obj.RHSE);
            obj.SimulationResult = obj.NonhydrostaticSolver.SPNPX\rhs(:);
            Err = obj.ExactResult(:) - obj.SimulationResult;
        end
        
        
    end
    
    methods ( Access = protected )
        
        function matGetFunction(obj)
            syms x y;
            obj.UE = sin(2*pi*x/20 + pi/2);
            obj.RHSE = diff(diff(obj.UE,x),x);
        end
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )           
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                % bottom elevation
                fphys2d{m}(:, :, 4) = -obj.H0;                
                %water depth
                fphys2d{m}(:,:,1) =  obj.mesh2d(m).x.*(20-obj.mesh2d(m).x).*obj.mesh2d(m).y.*(20-obj.mesh2d(m).y) - fphys2d{m}(:, :, 4);
%                 fphys2d{m}(:,:,1) =  sin(pi*obj.mesh2d(m).x/10).*sin(pi*obj.mesh2d(m).y/10) - fphys2d{m}(:, :, 4);
            end
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 20;
            outputIntervalNum = 1500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;                  
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK222;
%             option('EddyViscosityType') = enumEddyViscosity.Constant;
%             option('GOTMSetupFile') = obj.GotmFile;
%             option('equationType') = enumDiscreteEquation.Strong;
%             option('integralType') = enumDiscreteIntegral.GaussQuadrature;
%             option('outputType') = enumOutputFile.VTK;
%             option('ConstantEddyViscosityValue') = 0;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 1; 
            
        end
        
    end
end

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad, ...
    enumBoundaryCondition.ZeroGrad ];

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, ceil(obj.ChWidth/(obj.ChLength/M)), bctype);

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

