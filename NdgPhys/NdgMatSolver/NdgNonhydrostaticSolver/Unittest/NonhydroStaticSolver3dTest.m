classdef NonhydroStaticSolver3dTest < SWEBarotropic3d
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
    end
    
    properties
        U3d
        
        V3d
        
        W3d
        
        H
        
        Eta
        
        Sigma
        
        PSPX
        
        PSPY
      
        SPSPX
       
        SPSPY 
      
        SQPSPX
      
        SQPSPY
      
        PUPX
        
        PUPS
        
        PVPY
        
        PVPS
        
        PWPS        
        
    end
    
    properties
        
      ExactPSPX
      
      ExactPSPY
      
      ExactSPSPX
      
      ExactSPSPY 
      
      ExactSQPSPX
      
      ExactSQPSPY
      
      ExactPUPX
        
      ExactPUPS
        
      ExactPVPY
        
      ExactPVPS
        
      ExactPWPS
    end  
    
    properties
        Lambda = 20;
        A = 0.5;
    end
    
    methods
        function obj = NonhydroStaticSolver3dTest( N, Nz, M, Mz )
            % setup mesh domain
            [ mesh2d, mesh3d ] = makeChannelMesh( obj, N, Nz, M, Mz );
            obj.outputFieldOrder2d = [ 1 2 3 ];
            obj.outputFieldOrder3d = [1 2 3 4 5 6 7 8 9 10 11];
            obj.Nfield = 11;
            obj.Nvar = 11;
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, mesh3d );
            
            obj.matGetExactFunction;
            
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( mesh2d, mesh3d );
            
            obj.matGetExactSolution;
            
            obj.NonhydrostaticSolver.TestPartialDerivativeCalculation(obj, obj.fphys, obj.fphys2d);
            
        end
        
                
        function matGetExactSolution( obj )
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
            obj.ExactPSPX = eval(obj.PSPX);
            obj.ExactPSPY = eval(obj.PSPY);
            obj.ExactSPSPX = eval(obj.SPSPX);
            obj.ExactSPSPY = eval(obj.SPSPY);
            obj.ExactSQPSPX = eval(obj.SQPSPX);
            obj.ExactSQPSPY = eval(obj.SQPSPY);
            obj.ExactPUPX = eval(obj.PUPX);
            obj.ExactPUPS = eval(obj.PUPS);
            obj.ExactPVPY = eval(obj.PVPY);
            obj.ExactPVPS = eval(obj.PVPS);
            obj.ExactPWPS = eval(obj.PWPS);
        end
        
        
    end
    
    methods ( Access = protected )
        
        function matGetExactFunction( obj )
            syms x y z;
            obj.U3d = sin(2*pi*x/obj.Lambda)*sin(2*pi*y/obj.Lambda)*sin(2*pi*z);
            obj.V3d = sin(2*pi*x/obj.Lambda)*sin(2*pi*y/obj.Lambda)*sin(2*pi*z);
            obj.W3d = sin(2*pi*x/obj.Lambda)*sin(2*pi*y/obj.Lambda)*sin(2*pi*z);
            obj.H = obj.A * sin(2*pi*x/obj.Lambda) * sin(2*pi*y/obj.Lambda) + obj.H0;
            obj.Eta = obj.A * sin(2*pi*x/obj.Lambda) * sin(2*pi*y/obj.Lambda);
            obj.Sigma = (z - obj.Eta)/obj.H;
            obj.PSPX = diff(obj.Sigma, x);
            obj.PSPY = diff(obj.Sigma, y);
            obj.SPSPX = diff(diff(obj.Sigma, x), x);
            obj.SPSPY = diff(diff(obj.Sigma, y), y);
            obj.SQPSPX = obj.PSPX * obj.PSPX;
            obj.SQPSPY = obj.PSPY * obj.PSPY;
            obj.PUPX = diff(obj.U3d, x);
            obj.PUPS = diff(obj.U3d, z);
            obj.PVPY = diff(obj.V3d, y);
            obj.PVPS = diff(obj.V3d, z);
            obj.PWPS = diff(obj.W3d, z); 
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
                x = obj.mesh2d(m).x;
                y = obj.mesh2d(m).y;
                fphys2d{m}(:,:,1) =  eval(obj.H);
                x = obj.meshUnion(m).x;
                y = obj.meshUnion(m).y;
                z = obj.meshUnion(m).z;
                fphys{m}(:,:,1) = eval(obj.H).*eval(obj.U3d);
                fphys{m}(:,:,2) = eval(obj.H).*eval(obj.V3d);
                fphys{m}(:,:,11) = eval(obj.H).*eval(obj.W3d);
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
            option('outputType') = enumOutputFile.NetCDF;
            option('ConstantVerticalEddyViscosityValue') = 0.03;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 1; 
            
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
    [ 0, obj.ChLength ], [0, obj.ChWidth], M, ceil(obj.ChWidth/(obj.ChLength/M)), bctype);

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

