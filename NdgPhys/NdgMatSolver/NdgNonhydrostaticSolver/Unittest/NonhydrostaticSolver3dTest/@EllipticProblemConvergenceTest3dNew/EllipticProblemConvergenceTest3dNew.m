classdef EllipticProblemConvergenceTest3dNew < SWEBarotropic3d
    %ELLIPTICPROBLEMCONVERGENCETEST3DNEW 此处显示有关此类的摘要
    %   此处显示详细说明
    % Note: This case is exclusive for the test of convergence property
    
    properties
        ChLength = 2
        ChWidth = 2
        K11
        K22
        K13
        K23
        K33
        StiffMatrix
        SurfaceBoundaryEdgeType = 'Dirichlet'
    end
    
    properties (Constant)
        hcrit = 0.01
    end
    
    properties
        D
        zb
        Depth = 0.2
        Cexact
        symK11
        symK22
        symK13
        symK23
        symK33
        SecondDiffCexact
        ExactSolution
        SimulatedSolution
        RHS
    end
    
    methods
        function obj = EllipticProblemConvergenceTest3dNew(N, Nz, M, Mz)
            [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.varFieldIndex = [1 2 11];
            obj.K11 = ones(mesh3d.cell.Np, mesh3d.K);
            obj.K22 = obj.K11;
            obj.K13 = zeros(mesh3d.cell.Np, mesh3d.K);
            obj.K23 = zeros(mesh3d.cell.Np, mesh3d.K);
            obj.K33 = obj.K13.^2 + obj.K23.^2 + 1/obj.Depth/obj.Depth;
            obj.matGetFunction;
            obj.RHS = obj.matAssembleRightHandSide;
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.meshUnion );
        end
        
        function EllipticProblemSolve(obj)
            mesh3d = obj.meshUnion;
            x = mesh3d.x;
            y = mesh3d.y;
            z = mesh3d.z;
            obj.matGetFunction;
            obj.K11 = eval(obj.symK11) * ones(size(x));
            obj.K22 = obj.K11;
            obj.K13 = eval(obj.symK13) * ones(size(x));
            obj.K23 = eval(obj.symK23) * ones(size(x));
            obj.K33 = eval(obj.symK33) * ones(size(x));
            obj.RHS = obj.matAssembleRightHandSide;
            obj.AssembleGlobalStiffMatrix;
            
            NonhydroStiffMatrix = obj.NonhydrostaticSolver.TestAssembleGlobalStiffMatrix( obj, obj.fphys, obj.fphys2d );
            
            Matrix = full(obj.StiffMatrix - NonhydroStiffMatrix);
            
            disp("============For Matrix================")
            disp("The maximum difference is:")
            disp(max(max(Matrix)))
            disp("The minimum difference is:")
            disp(min(min(Matrix)))
            disp("============End Matrix================")
            
            obj.SimulatedSolution = NonhydroStiffMatrix\(obj.RHS(:));
            
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
        
        AssembleGlobalStiffMatrix(obj);
        
    end
    methods(Access = protected)
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            syms x y;
            %variable water depth
            obj.zb = 0*x + 0*y + 0.2;
            obj.D = 0*x + 0*y + obj.Depth;
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh2d = obj.meshUnion(m).mesh2d;
                mesh3d = obj.meshUnion(m);
                x = mesh2d.x;
                y = mesh2d.y;
                fphys2d{m} = zeros( mesh2d.cell.Np, mesh2d.K, obj.Nfield2d );
                fphys{m} = zeros( mesh3d.cell.Np, mesh3d.K, obj.Nfield );
                fphys2d{m}(:, :, 4) = eval(obj.zb);
                %water depth
                fphys2d{m}(:,:,1) = eval(obj.D);
                %                  fphys2d{m}(:,:,1) = 2.89677;
            end
        end
        
        function rhs = matAssembleRightHandSide( obj )
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
            Temprhs = eval(obj.SecondDiffCexact);
            for i = 1:size(Temprhs,2)
                rhs(:,i) = diag(obj.meshUnion.J(:,1))*obj.meshUnion.cell.M * Temprhs(:,i);
            end
        end
        
        function matGetFunction(obj)
            syms x y z nx ny nz;
            obj.Cexact = sin(-pi/2*x) * sin(-pi/2*y) * sin(-pi/2*z);
            obj.symK11 = 0*x + 0*y + 1;
            obj.symK22 = 0*x + 0*y + 1;
            obj.symK13 =1*( -1/obj.D*(diff(obj.D,x)) - z/obj.D*diff(obj.D,x));
            obj.symK23 =1*(-1/obj.D*(diff(obj.D,y)) - z/obj.D*diff(obj.D,y));
            obj.symK33 = obj.symK13^2 + obj.symK23^2 + 1./obj.D./obj.D;
            
            obj.SecondDiffCexact = diff(diff(obj.Cexact,x) + obj.symK13*diff(obj.Cexact, z),x) + ...
                diff(diff(obj.Cexact,y) + obj.symK23*diff(obj.Cexact, z),y) + ...
                diff(obj.symK13*diff(obj.Cexact,x) + obj.symK23*diff(obj.Cexact,y) + obj.symK33*diff(obj.Cexact, z),z);
            
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
            obj.ExactSolution = eval(obj.Cexact);
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
            outputIntervalNum = 7500;
            option('startTime') = 0.0;
            option('finalTime') = 1.0;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = 1.0/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 20;
            option('temporalDiscreteType') = enumTemporalDiscrete.SSPRK22;
            option('VerticalEddyViscosityType') = enumSWEVerticalEddyViscosity.None;
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.VTK;
            option('limiterType') = enumLimiter.Vert;
            option('ConstantVerticalEddyViscosityValue') = 0.0001;
            option('HorizontalEddyViscosityType') = enumSWEHorizontalEddyViscosity.None;
            option('ConstantHorizontalEddyViscosityValue') = 1.0;
            option('BottomBoundaryEdgeType') = enumBottomBoundaryEdgeType.Neumann;
            option('CoriolisType') = enumSWECoriolis.None;
        end
    end
end

