classdef EllipticProblemInHorizontalDirectionTest3d < SWEBarotropic3d
    
    %> Note: This solver is used for test purpose of operator
    %> $\frac{\partial^2p}{\partial \x^2} + \frac{\partial^2p}{\partial \y^2}$, the analytical
    %> solution $p=\pi^5 sin\left (-\frac{\pi}{10}x\right ) + \pi sin\left (-\frac{\pi}{10}y\right )$ is used 
    properties
        ChLength = 20
        ChWidth = 0.4
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
    end
    
    properties(Constant)
        hcrit = 1
    end
    
    
    methods
        function obj = EllipticProblemInHorizontalDirectionTest3d(N, Nz, M, Mz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.meshUnion = obj.mesh3d;
            obj.Nmesh = 1;
            obj.matGetFunction;
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.mesh3d );
            obj.RHS = obj.matAssembleRightHandSide;
            obj.matAssembleGlobalStiffMatrix;
        end
        function EllipticProblemSolve(obj)
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            obj.ExactSolution = eval(obj.Cexact);
            obj.SimulatedSolution = obj.StiffMatrix\obj.RHS(:);
            disp("The maximum difference is:");
            disp(max(max(obj.ExactSolution(:)-obj.SimulatedSolution)));
            disp("The minimum difference is:");
            disp(min(min(obj.ExactSolution(:)-obj.SimulatedSolution)));
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
        
        function matGetFunction(obj)
            syms x y;
            obj.Cexact = pi^5*sin(-pi/10*x) + pi * sin(-pi/10*y) + 30;
            obj.SecondDiffCexact = diff(diff(obj.Cexact, x), x) + diff(diff(obj.Cexact, y),y);
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