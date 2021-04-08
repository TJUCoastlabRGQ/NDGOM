classdef EllipticProblemInHorizontalDirectionTest3d < SWEBarotropic3d
    
    %> Note: This solver is used for test purpose of operator
    %> $\frac{\partial^2p}{\partial \x^2} + \frac{\partial^2p}{\partial \y^2}$, the analytical
    %> solution $p=\pi^5 sin\left (-\frac{\pi}{10}x\right ) + \pi sin\left (-\frac{\pi}{10}y\right )$ is used
    properties
        ChLength = 20
        ChWidth = 4
    end
    
    properties
        ExactSolution
        
        ExactRHS
        
        StiffMatrix
        
        RHS
        
        SimulatedSolution
        
        Cexact
        
        DiffCexact
        
        SecondDiffCexact
        
        DirichExact
        
        NewmannExact
        
        DirichletData
    end
    
    properties(Constant)
        hcrit = 1
        D11 = 1
        D22 = 1
    end
    
    
    methods
        function obj = EllipticProblemInHorizontalDirectionTest3d(N, Nz, M, Mz)
            [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.matGetFunction;
            
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.meshUnion );
            obj.RHS = obj.matAssembleRightHandSide;
            
            x = obj.NonhydrostaticSolver.BoundaryEdge.xb;
            y = obj.NonhydrostaticSolver.BoundaryEdge.yb;
            z = obj.NonhydrostaticSolver.BoundaryEdge.zb;
            obj.DirichletData = eval(obj.Cexact);
            
            obj.matAssembleGlobalStiffMatrix;
            
        end
        function EllipticProblemSolve(obj)
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
            obj.ExactSolution = eval(obj.Cexact);
            obj.SimulatedSolution = obj.StiffMatrix\obj.RHS(:);
%             disp("The condition number is:");
%             disp(condest(obj.StiffMatrix));
%             disp("The maximum difference is:");
%             disp(max(max(obj.ExactSolution(:)-obj.SimulatedSolution)));
%             disp("The minimum difference is:");
%             disp(min(min(obj.ExactSolution(:)-obj.SimulatedSolution)));
        end
        
    end
    
    methods( Hidden )
        function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys )
            %an empty function
        end
    end
    
    methods(Access = protected)
        
        matAssembleGlobalStiffMatrix(obj);
        
        rhs = matAssembleRightHandSide(obj);
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
            end
        end
        
        function matGetFunction(obj)
            syms x y z;
            obj.Cexact = sin(-pi/10*x)*sin(-pi/10*y);
            obj.SecondDiffCexact = diff(diff(obj.Cexact, x), x) + diff(diff(obj.Cexact, y),y);
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
    [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

end