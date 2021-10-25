classdef EllipticProblemConvergenceTest3dNew < EllipticProblemMatrixAssembleTest3dNew
    %ELLIPTICPROBLEMCONVERGENCETEST3DNEW 此处显示有关此类的摘要
    %   此处显示详细说明
    
    methods
        function obj = EllipticProblemConvergenceTest3dNew(N, Nz, M, Mz)
            obj = obj@EllipticProblemMatrixAssembleTest3dNew(N, Nz, M, Mz);
        end
        
        function EllipticProblemSolve(obj)
            mesh3d = obj.meshUnion;
            obj.K11 = ones(mesh3d.cell.Np, mesh3d.K);
            obj.K22 = obj.K11;
            obj.K13 = 0.04 * ones(mesh3d.cell.Np, mesh3d.K);
            obj.K23 = 0.04 * ones(mesh3d.cell.Np, mesh3d.K);
            obj.K33 = (obj.K13.^2 + obj.K23.^2 + 1/obj.Depth/obj.Depth);
            obj.matGetFunction;
            obj.RHS = obj.matAssembleRightHandSide;
            obj.NonhydrostaticSolver.TestNewFormGlobalStiffMatrix( obj, obj.fphys );
            obj.matAssembleGlobalStiffMatrixWithBCsImposed;
            obj.matAssembleGlobalStiffMatrixWithSurfaceBCsImposed;
            obj.matAssembleGlobalStiffMatrixWithBottomBCsImposed;
            
            obj.SimulatedSolution = obj.NonhydrostaticSolver.GlobalStiffMatrix\obj.RHS(:);
            
            disp("============For solution================")
            disp("The maximum value is:")
            disp(max(max(obj.ExactSolution)))
            disp("The minimum value is:")
            disp(min(min(obj.ExactSolution)))
            disp("The maximum difference is:")
            disp(max(max(obj.ExactSolution(:) - obj.SimulatedSolution(:))));
            disp("The minimum difference is:")
            disp(min(min(obj.ExactSolution(:) - obj.SimulatedSolution(:))));
            disp("============End solution================")
        end
        
    end
    methods(Access = protected)
        
        function matGetFunction(obj)
            syms x y z nx ny nz;
            
            obj.Cexact = sin(-pi/2*z)+sin(pi/2*x)+sin(pi/2*y);    

            obj.SecondDiffCexact = diff(obj.K11(1)*diff(obj.Cexact,x) + obj.K13(1)*diff(obj.Cexact, z),x) + ...
                diff(obj.K22(1)*diff(obj.Cexact,y) + obj.K23(1)*diff(obj.Cexact, z),y) + ...
                diff(obj.K13(1)*diff(obj.Cexact,x) + obj.K23(1)*diff(obj.Cexact,y) + obj.K33(1)*diff(obj.Cexact, z),z);
            
            obj.NewmannCexact = (obj.K11(1)*diff(obj.Cexact,x) + obj.K13(1)*diff(obj.Cexact, z))*nx + ...
                (obj.K22(1)*diff(obj.Cexact,y) + obj.K23(1)*diff(obj.Cexact, z))*ny + ...
                (obj.K13(1)*diff(obj.Cexact,x) + obj.K23(1)*diff(obj.Cexact,y) + obj.K33(1)*diff(obj.Cexact, z))*nz;
            
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
            obj.ExactSolution = eval(obj.Cexact);
            
            x = obj.meshUnion.BoundaryEdge.xb;
            y = obj.meshUnion.BoundaryEdge.yb;
            z = obj.meshUnion.BoundaryEdge.zb;
            nx = obj.meshUnion.BoundaryEdge.nx;
            ny = obj.meshUnion.BoundaryEdge.ny;
            nz = obj.meshUnion.BoundaryEdge.nz;
            obj.DirichletData = eval(obj.Cexact);
            obj.NewmannData = eval(obj.NewmannCexact);
            
            x = obj.meshUnion.mesh2d.x;
            y = obj.meshUnion.mesh2d.y;
            z = max(max(obj.meshUnion.z))*ones(size(x));
            nx = obj.meshUnion.SurfaceBoundaryEdge.nx;
            ny = obj.meshUnion.SurfaceBoundaryEdge.ny;
            nz = obj.meshUnion.SurfaceBoundaryEdge.nz;
            obj.SurfaceDirichletData = eval(obj.Cexact);
            obj.SurfaceNewmannData = eval(obj.NewmannCexact);
            
            x = obj.meshUnion.mesh2d.x;
            y = obj.meshUnion.mesh2d.y;
            z = min(min(obj.meshUnion.z))*ones(size(x));
            nx = obj.meshUnion.BottomBoundaryEdge.nx;
            ny = obj.meshUnion.BottomBoundaryEdge.ny;
            nz = obj.meshUnion.BottomBoundaryEdge.nz;
            obj.BottomDirichletData = eval(obj.Cexact);
            obj.BottomNewmannData = eval(obj.NewmannCexact);
        end
        
        function [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz)
            bctype = [ ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
                enumBoundaryCondition.SlipWall, ...
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

