classdef EllipticMixedParticalDerivativeTest3d < SWEBarotropic3d
    %> For term $\frac{\partial }{\partial x}\left ( \frac{\partial p}{\partial y}\right )\right ) = f$
    
    properties
        ChLength = 2
        ChWidth = 1
    end
    
    properties
        ExactSolution
        
        ExactRHS
        
        StiffMatrix
        
        RHS
        
        SimulatedSolution
        
        Cexact
        
        DiffCexact
        
        MixedSecondDiffCexact
        
        DirichExact
        
        NewmannExact
                
        SurfaceDirichletData
        
        BottomBoundaryDirichletData
        
        BoundaryDirichletData
    end
    
    properties(Constant)
        hcrit = 1
    end
    
    methods
        
        function obj = EllipticMixedParticalDerivativeTest3d(N, Nz, M, Mz)
            [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.matGetFunction;
            
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.meshUnion );
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
            obj.RHS = eval(obj.MixedSecondDiffCexact);
            
            x = obj.NonhydrostaticSolver.BoundaryEdge.xb;
            y = obj.NonhydrostaticSolver.BoundaryEdge.yb;
            z = obj.NonhydrostaticSolver.BoundaryEdge.zb;
            obj.BoundaryDirichletData = eval(obj.Cexact);
            
            x = obj.meshUnion.mesh2d.x;
            y = obj.meshUnion.mesh2d.y;
            z = obj.meshUnion.mesh2d.z;
            obj.SurfaceDirichletData = eval(obj.Cexact);
            
            x = obj.meshUnion.mesh2d.x;
            y = obj.meshUnion.mesh2d.y;
            z = -1 * ones(size(obj.meshUnion.mesh2d.z));
            obj.BottomBoundaryDirichletData = eval(obj.Cexact);            
            
            
            obj.AssembleGlobalStiffMatrix;
            
        end
        function EllipticProblemSolve(obj)
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            z = obj.meshUnion.z;
            obj.ExactSolution = eval(obj.Cexact);
            obj.SimulatedSolution = obj.StiffMatrix\obj.RHS(:);
            disp("The condition number is:");
            disp(condest(obj.StiffMatrix));
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
            obj.Cexact = sin(-pi*x)*sin(pi*y)*sin(-pi/2*z);
            obj.MixedSecondDiffCexact = diff(diff(obj.Cexact, z), x) + diff(diff(obj.Cexact, z),y);
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
[ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
[ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end