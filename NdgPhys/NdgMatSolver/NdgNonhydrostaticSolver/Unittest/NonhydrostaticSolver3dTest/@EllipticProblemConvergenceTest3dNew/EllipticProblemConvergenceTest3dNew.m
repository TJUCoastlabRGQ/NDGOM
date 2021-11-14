classdef EllipticProblemConvergenceTest3dNew < EllipticProblemMatrixAssembleTest3dNew
    %ELLIPTICPROBLEMCONVERGENCETEST3DNEW 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        ConstantCoe = 'False'
    end
    
    methods
        function obj = EllipticProblemConvergenceTest3dNew(N, Nz, M, Mz)
            obj = obj@EllipticProblemMatrixAssembleTest3dNew(N, Nz, M, Mz);
        end
        
        function EllipticProblemSolve(obj)
            mesh3d = obj.meshUnion;
            x = mesh3d.x;
            y = mesh3d.y;
            z = mesh3d.z;
            obj.matGetFunction;
            obj.K11 = eval(obj.symK11);
            if numel(obj.K11) == 1
                obj.K11 = obj.K11 * ones(size(x));
            end
            obj.K22 = obj.K11;
            obj.K13 = eval(obj.symK13);
            if numel(obj.K13) == 1
                obj.ConstantCoe = 'True';
                obj.K13 = obj.K13 * ones(size(x));
            end
            obj.K23 = eval(obj.symK23);
            if numel(obj.K23) == 1
                obj.K23 = obj.K23 * ones(size(x));
            end            
            obj.K33 = eval(obj.symK33);
            if numel(obj.K33) == 1
                obj.K33 = obj.K33 * ones(size(x));
            end            
            obj.RHS = obj.matAssembleRightHandSide;
%             obj.NonhydrostaticSolver.TestNewFormGlobalStiffMatrix( obj, obj.fphys, obj.fphys2d );
            obj.NonhydrostaticSolver.AssembleStiffMatrix( obj, obj.fphys, obj.fphys2d );
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
        
                %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            syms x y;
            %variable water depth
            obj.D = 1*0.1*cos(pi*x)*cos(pi*y)+obj.Depth;
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                x = obj.mesh2d(m).x;
                y = obj.mesh2d(m).y;
                fphys2d{m}(:,:,1) = eval(obj.D);
            end
        end
        
        function matGetFunction(obj)
            syms x y z nx ny nz;
            obj.Cexact = sin(-pi/2*z)*sin(1/2*x)*sin(1/2*y);             
            obj.symK11 = 0*x + 0*y + 1;
            obj.symK22 = 0*x + 0*y + 1;
            obj.symK13 =1*( -1/obj.D*(diff(obj.D,x)) - z/obj.D*diff(obj.D,x)) + 0*0.04;
            obj.symK23 =1*(-1/obj.D*(diff(obj.D,y)) - z/obj.D*diff(obj.D,y)) + 0*0.04;
            obj.symK33 = obj.symK13^2 + obj.symK23^2 + 1./obj.D./obj.D;
            
            obj.SecondDiffCexact = diff(diff(obj.Cexact,x) + obj.symK13*diff(obj.Cexact, z),x) + ...
                diff(diff(obj.Cexact,y) + obj.symK23*diff(obj.Cexact, z),y) + ...
                diff(obj.symK13*diff(obj.Cexact,x) + obj.symK23*diff(obj.Cexact,y) + obj.symK33*diff(obj.Cexact, z),z) - ...
                diff(obj.Cexact,x)*diff(obj.symK13,z) - diff(obj.Cexact,z)*obj.symK13*diff(obj.symK13,z) - ...
                diff(obj.Cexact,y)*diff(obj.symK23,z) - diff(obj.Cexact,z)*obj.symK23*diff(obj.symK23,z);
            
            obj.NewmannCexact = (diff(obj.Cexact,x) + obj.symK13*diff(obj.Cexact, z))*nx + ...
                (diff(obj.Cexact,y) + obj.symK23*diff(obj.Cexact, z))*ny + ...
                (obj.symK13*diff(obj.Cexact,x) + obj.symK23*diff(obj.Cexact,y) + obj.symK33*diff(obj.Cexact, z))*nz;
            
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

