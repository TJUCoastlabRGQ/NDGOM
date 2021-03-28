classdef EllipticProblemInVerticalDirectionTest3d < SWEAbstract3d
    
    %> Note: This solver is used for test purpose of operator
    %> $\frac{\partial^2p}{\partial \sigma^2}$, the analytical
    %> solution $p=sin\left (-\frac{\pi}{2}z\right )$ is used 
    
    properties
        ChLength = 2
        ChWidth = 2
        Nfield2d = 1
        Nvar2d = 1
        varFieldIndex2d = 1
    end
    
    properties
        
        StiffMatrix
        
        RHS
        
        ExactSolution
        
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
    
    properties
        outputFieldOrder2d = []
        outputFieldOrder3d =  1
        Nfield = 1
        Nvar = 1
        varFieldIndex = 1
    end
    
    methods
        function obj = EllipticProblemInVerticalDirectionTest3d(N, Nz, M, Mz)
            [obj.mesh2d, obj.mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.meshUnion = obj.mesh3d;
            obj.Nmesh = 1;
            obj.matGetFunction;
            obj.DirichExact = zeros(obj.mesh2d.cell.Np,1);
            obj.NewmannExact = zeros(obj.mesh2d.cell.Np,1);
            [ obj.StiffMatrix, LRHS, RRHS ] = obj.matAssembleGlobalStiffMatrix;
            obj.RHS = obj.matAssembleRightHandSide;
            % The direction vector is considered in the initialization stage
            obj.RHS(:,1) = obj.RHS(:,1) - LRHS;
            % Minus or plus is considered when assemble RRHS
            obj.RHS(:,end) = obj.RHS(:,end) + RRHS;            
        end
        function EllipticProblemSolve(obj)
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
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
        
        [ Matrix, Lrhs, Rrhs ] = matAssembleGlobalStiffMatrix(obj);
        
        rhs = matAssembleRightHandSide(obj);
        
        function matGetFunction(obj)
            syms z;
            obj.Cexact = sin(-pi/2*z);
            obj.DiffCexact = diff(obj.Cexact, z);
            obj.SecondDiffCexact = diff(obj.DiffCexact, z);
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