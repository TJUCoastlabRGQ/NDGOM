classdef FirstOrderProblemInVerticalTest3d < SWEBarotropic3d
    
    %> Note: This solver is used for test purpose of operator
    %> $\frac{\partial p}{\partial \sigma }$, the analytical
    %> solution $p=sin\left (-\frac{\pi}{2}z\right )$ is used   
    
    properties
        ChLength = 3
        ChWidth = 3
    end
    
    properties
        
        StiffMatrix
        
        RHS
        
        ExactSolution
        
        SimulatedSolution
        
        Cexact
        
        DiffCexact
                
        DirichExact
        
        SurfaceBoundaryEdgeType = "Dirichlet"
        
    end
    
    properties(Constant)
        hcrit = 1
    end
    
    methods
        function obj = FirstOrderProblemInVerticalTest3d(N, Nz, M, Mz)
            [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.varFieldIndex = [1 2 11];
            obj.matGetFunction;
            z = zeros(obj.mesh2d.cell.Np,1);
            obj.DirichExact = eval(obj.Cexact);
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.meshUnion );
            [ obj.StiffMatrix, LRHS ] = obj.matAssembleGlobalStiffMatrix;
            obj.RHS = obj.matAssembleRightHandSide;
            % The direction vector is considered in the initialization stage
            obj.RHS(:,1) = obj.RHS(:,1) + LRHS;
        end
        function EllipticProblemSolve(obj)
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
            obj.ExactSolution = eval(obj.Cexact);
%             obj.SimulatedSolution = obj.StiffMatrix\obj.RHS(:);
            obj.SimulatedSolution = obj.NonhydrostaticSolver.PNPS\obj.RHS(:);
            disp("The condition number is:")
            disp(condest(sparse(obj.StiffMatrix)));
            disp("============For solution===================");
            disp("The maximum difference is:");
            disp(max(max(obj.ExactSolution(:)-obj.SimulatedSolution)));
            disp("The minimum difference is:");
            disp(min(min(obj.ExactSolution(:)-obj.SimulatedSolution)));
            disp("============End solution===================");
            disp("============For stiff matrix===============");
            disp("The maximum difference is:");
            disp(max(max(obj.StiffMatrix - obj.NonhydrostaticSolver.PNPS)));
            disp("The minimum difference is:");
            disp(min(min(obj.StiffMatrix - obj.NonhydrostaticSolver.PNPS)));
            disp("============End stiff matrix==============")
            
        end
        
    end
    
    methods( Hidden )
        function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys )
            %an empty function
        end
    end
    
    methods(Access = protected)
        
        [ Matrix, Lrhs ] = matAssembleGlobalStiffMatrix(obj);
        
        rhs = matAssembleRightHandSide(obj);
        
        function matGetFunction(obj)
            syms z;
            obj.Cexact = sin(-pi/2*z);
            obj.DiffCexact = diff(obj.Cexact, z);
        end
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
            end
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

function [mesh2d, mesh3d] = makeChannelMesh( obj, N, Nz, M, Mz )

bctype = [ ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall, ...
    enumBoundaryCondition.SlipWall ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], 1, 1, bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end