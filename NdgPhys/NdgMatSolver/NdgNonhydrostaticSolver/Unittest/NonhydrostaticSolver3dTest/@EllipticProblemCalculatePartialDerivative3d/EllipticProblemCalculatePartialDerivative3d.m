classdef EllipticProblemCalculatePartialDerivative3d < SWEBarotropic3d
    %ELLIPTICPROBLEMCALCULATEPARTIALDERIVATIVE3D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        matPSPX
        matPSPY
        matSQPSPX
        matSQPSPY
        matPUPX
        matPUPY
        matPVPX
        matPVPY
        matPUPS
        matPVPS
        matPWPS
        matWnew
        matUnew
        matVnew
        matPHPX
        matPHPY
    end
    
    properties
        ChLength = 10
        ChWidth = 10
        SurfaceBoundaryEdgeType = 'Dirichlet'
    end
    
    properties(Constant)
        hcrit = -0.01
    end
    
    methods
        function obj = EllipticProblemCalculatePartialDerivative3d(N, Nz, M, Mz)
            [mesh2d, mesh3d] = obj.makeChannelMesh(N, Nz, M, Mz);
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.varFieldIndex = [1 2 11];
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.meshUnion );
        end
        
        function EllipticProblemSolve(obj)
            
           mesh2d = obj.meshUnion(1).mesh2d;
           
           edge3d = obj.meshUnion(1).BoundaryEdge;
           % For Hu
           obj.fext3d{1}(:,:,1) = rand(edge3d.Nfp, edge3d.Ne);
           % For Hv
           obj.fext3d{1}(:,:,2) = rand(edge3d.Nfp, edge3d.Ne);
           % For H
           obj.fext3d{1}(:,:,3) = rand(edge3d.Nfp, edge3d.Ne);
           % For Hw
           obj.fext3d{1}(:,:,11) = rand(edge3d.Nfp, edge3d.Ne);
           
           obj.frhs2d{1} = rand(mesh2d.cell.Np, mesh2d.K);
           
           obj.NonhydrostaticSolver.TestPartialDerivativeCalculation(obj, obj.fphys, obj.fphys2d );
           
           obj.matCalculatePartialDerivative;                     
           
           disp('===============For PSPX=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPSPX)))
           disp('The minimum value is')
           disp(min(min(obj.matPSPX)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPSPX - obj.NonhydrostaticSolver.PSPX)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPSPX - obj.NonhydrostaticSolver.PSPX)))
           disp('===============End PSPX=================')
           
           disp('===============For PSPY=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPSPY)))
           disp('The minimum value is')
           disp(min(min(obj.matPSPY)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPSPY - obj.NonhydrostaticSolver.PSPY)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPSPY - obj.NonhydrostaticSolver.PSPY)))
           disp('===============End PSPY=================')
           
           disp('===============For SQPSPX=================')
           disp('The maximum value is:')
           disp(max(max(obj.matSQPSPX)))
           disp('The minimum value is')
           disp(min(min(obj.matSQPSPX)))
           disp('The maximum difference is:')
           disp(max(max(obj.matSQPSPX - obj.NonhydrostaticSolver.SQPSPX)))
           disp('The minimum difference is:')
           disp(min(min(obj.matSQPSPX - obj.NonhydrostaticSolver.SQPSPX)))
           disp('===============End SQPSPX=================')
           
           disp('===============For SQPSPY=================')
           disp('The maximum value is:')
           disp(max(max(obj.matSQPSPY)))
           disp('The minimum value is')
           disp(min(min(obj.matSQPSPY)))
           disp('The maximum difference is:')
           disp(max(max(obj.matSQPSPY - obj.NonhydrostaticSolver.SQPSPY)))
           disp('The minimum difference is:')
           disp(min(min(obj.matSQPSPY - obj.NonhydrostaticSolver.SQPSPY)))
           disp('===============End SQPSPY=================')
           
           disp('===============For PUPX=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPUPX)))
           disp('The minimum value is')
           disp(min(min(obj.matPUPX)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPUPX - obj.NonhydrostaticSolver.PUPX)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPUPX - obj.NonhydrostaticSolver.PUPX)))
           disp('===============End PUPX=================')
           
           disp('===============For PUPY=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPUPY)))
           disp('The minimum value is')
           disp(min(min(obj.matPUPY)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPUPY - obj.NonhydrostaticSolver.PUPY)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPUPY - obj.NonhydrostaticSolver.PUPY)))
           disp('===============End PUPY=================')  
           
           disp('===============For PVPX=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPVPX)))
           disp('The minimum value is')
           disp(min(min(obj.matPVPX)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPVPX - obj.NonhydrostaticSolver.PVPX)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPVPX - obj.NonhydrostaticSolver.PVPX)))
           disp('===============End PVPX=================')   
           
           disp('===============For PVPY=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPVPY)))
           disp('The minimum value is')
           disp(min(min(obj.matPVPY)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPVPY - obj.NonhydrostaticSolver.PVPY)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPVPY - obj.NonhydrostaticSolver.PVPY)))
           disp('===============End PVPY=================') 
           
           disp('===============For PUPS=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPUPS)))
           disp('The minimum value is')
           disp(min(min(obj.matPUPS)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPUPS - obj.NonhydrostaticSolver.PUPS)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPUPS - obj.NonhydrostaticSolver.PUPS)))
           disp('===============End PUPS=================')
           
           disp('===============For PVPS=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPVPS)))
           disp('The minimum value is')
           disp(min(min(obj.matPVPS)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPVPS - obj.NonhydrostaticSolver.PVPS)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPVPS - obj.NonhydrostaticSolver.PVPS)))
           disp('===============End PVPS=================')  
           
           disp('===============For PWPS=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPWPS)))
           disp('The minimum value is')
           disp(min(min(obj.matPWPS)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPWPS - obj.NonhydrostaticSolver.PWPS)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPWPS - obj.NonhydrostaticSolver.PWPS)))
           disp('===============End PWPS=================')  
           
           disp('===============For Wnew=================')
           disp('The maximum value is:')
           disp(max(max(obj.matWnew)))
           disp('The minimum value is')
           disp(min(min(obj.matWnew)))
           disp('The maximum difference is:')
           disp(max(max(obj.matWnew - obj.NonhydrostaticSolver.Wnew)))
           disp('The minimum difference is:')
           disp(min(min(obj.matWnew - obj.NonhydrostaticSolver.Wnew)))
           disp('===============End Wnew=================')   
           
           disp('===============For Unew=================')
           disp('The maximum value is:')
           disp(max(max(obj.matUnew)))
           disp('The minimum value is')
           disp(min(min(obj.matUnew)))
           disp('The maximum difference is:')
           disp(max(max(obj.matUnew - obj.NonhydrostaticSolver.Unew)))
           disp('The minimum difference is:')
           disp(min(min(obj.matUnew - obj.NonhydrostaticSolver.Unew)))
           disp('===============End Unew=================')
           
           disp('===============For Vnew=================')
           disp('The maximum value is:')
           disp(max(max(obj.matVnew)))
           disp('The minimum value is')
           disp(min(min(obj.matVnew)))
           disp('The maximum difference is:')
           disp(max(max(obj.matVnew - obj.NonhydrostaticSolver.Vnew)))
           disp('The minimum difference is:')
           disp(min(min(obj.matVnew - obj.NonhydrostaticSolver.Vnew)))
           disp('===============End Vnew=================')   
           
           disp('===============For PHPX=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPHPX)))
           disp('The minimum value is')
           disp(min(min(obj.matPHPX)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPHPX - obj.NonhydrostaticSolver.PHPX)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPHPX - obj.NonhydrostaticSolver.PHPX)))
           disp('===============End PHPX=================')  
           
           disp('===============For PHPY=================')
           disp('The maximum value is:')
           disp(max(max(obj.matPHPY)))
           disp('The minimum value is')
           disp(min(min(obj.matPHPY)))
           disp('The maximum difference is:')
           disp(max(max(obj.matPHPY - obj.NonhydrostaticSolver.PHPY)))
           disp('The minimum difference is:')
           disp(min(min(obj.matPHPY - obj.NonhydrostaticSolver.PHPY)))
           disp('===============End PHPY=================')             
        end
    end
    methods(Access = protected)
                %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                % Hu
                fphys{m}(:,:,1) = 10 * rand(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K);
                % Hv
                fphys{m}(:,:,2) = 3 * rand(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K);
                % Hw
                fphys{m}(:,:,11) = rand(obj.meshUnion(m).cell.Np, obj.meshUnion(m).K) + ones( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K );
                % H
                fphys2d{m}(:,:,1) = 5 * rand( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K ) + ones( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K );
                % Z
                fphys2d{m}(:,:,4) = -rand(1)*obj.mesh2d(1).x - rand(1)*obj.mesh2d(1).y;
            end
        end
        
        function [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz)
            
            bctype = [ ...
                enumBoundaryCondition.ClampedDepth, ...
                enumBoundaryCondition.Clamped, ...
                enumBoundaryCondition.ClampedVel, ...
                enumBoundaryCondition.SlipWall ];
            mesh2d = makeUniformQuadMesh( N, ...
                [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], M, M, bctype);
            cell = StdPrismQuad( N, Nz );
            zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
            mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
            mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
            mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
            mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
            mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
            mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 1.75;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 5;
            option('AdvDiffVerticalDiffusionType') = enumVerticalDiffusion.None;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('AdvDiffHorizontalDiffusionType') = enumHorizontalDiffusion.Constant;
            option('AdvDiffConstantHorizontalDiffusionValue') = 1;
        end
    end
end

