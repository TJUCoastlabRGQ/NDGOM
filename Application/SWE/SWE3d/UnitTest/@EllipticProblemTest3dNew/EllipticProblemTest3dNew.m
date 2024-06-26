classdef EllipticProblemTest3dNew < SWEBarotropic3d
    
    %> Note: This solver is used for test purpose of operator
    %> $$\nabla\cdot\mathbf{K}\nabla p=f$$. Here, $\mathbf {K}$ is
    %> positive definite, and is given by,
    %> $$\mathbf K = \left [\begin{matrix}1&&0&&k_{13}\\0&&1&&k_{23}\\
    %> k_{31}&&k_{32}&&k_{33}\end {matrix}\right]$$. $k_{33}$ is calculated
    %> as $k_{13}^2+k_{23}^2+\frac{1}{D^2}$
    
    properties
        ChLength = 1
        ChWidth = 1
        Depth = 1
        TempK11 = 1
        TempK22 = 1
        TempK13 = 0.04
        TempK23 = 0.04
        TempK33
        K11
        K22
        K13
        K23
        K33
        SurfaceBoundaryEdgeType = 'Dirichlet'
        BottomBoundaryEdgeType = 'Newmann'        
    end
    
    properties
        
        StiffMatrix
        
        RHS
        
        ExactSolution
        
        SimulatedSolution
        
        Cexact
        
        SecondDiffCexact
        
        NewmannCexact
        
        DirichletData
        
        SurfaceDirichletData
        
        SurfaceNewmannData
        
        BottomDirichletData
        
        BottomNewmannData
        
        NewmannData
    end
    
    properties(Constant)
        hcrit = 0.01
    end
    
    
    methods
        function obj = EllipticProblemTest3dNew(N, Nz, M, Mz)
            [mesh2d, mesh3d] = makeChannelMesh(obj, N, Nz, M, Mz);
            obj.initPhysFromOptions( mesh2d, mesh3d );
            obj.K11 = obj.TempK11 * ones(obj.meshUnion.cell.Np, obj.meshUnion.K);
            obj.K22 = obj.TempK22 * ones(obj.meshUnion.cell.Np, obj.meshUnion.K);
            obj.K13 = obj.TempK13 * ones(obj.meshUnion.cell.Np, obj.meshUnion.K);
            obj.K23 = obj.TempK23 * ones(obj.meshUnion.cell.Np, obj.meshUnion.K);
            obj.K33 = obj.K13.^2 + obj.K23.^2 + 1/obj.Depth/obj.Depth;
%             obj.K33 = zeros(obj.meshUnion.cell.Np, obj.meshUnion.K);
            obj.TempK33 = obj.K33(1);
            obj.matGetFunction;
            obj.NonhydrostaticSolver = NdgQuadratureFreeNonhydrostaticSolver3d( obj, obj.meshUnion );
            obj.RHS = obj.matAssembleRightHandSide;
            obj.AssembleGlobalStiffMatrix;
        end
        function EllipticProblemSolve(obj)
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
            
            [ ~ ]= obj.NonhydrostaticSolver.TestNewFormGlobalStiffMatrix( obj, obj.fphys );       
            
            warning('off');
            [ obj.RHS, obj.NonhydrostaticSolver.GlobalStiffMatrix ] = mxAssembleGlobalStiffMatrixWithBCsImposed(...
                obj.NonhydrostaticSolver.GlobalStiffMatrix, obj.RHS, obj.DirichletData,...
                struct(obj.meshUnion.BoundaryEdge), struct(obj.meshUnion.cell), struct(obj.meshUnion), obj.K13, ...
                obj.K23, obj.K33, int8(obj.meshUnion.BoundaryEdge.ftype), obj.NewmannData);
            [ obj.RHS, ~ ] = mxAssembleGlobalStiffMatrixWithSurfaceBCsImposed(...
                obj.NonhydrostaticSolver.GlobalStiffMatrix, obj.RHS, obj.SurfaceDirichletData,...
                struct(obj.meshUnion.SurfaceBoundaryEdge), struct(obj.meshUnion.cell), struct(obj.meshUnion), obj.K13, ...
                obj.K23, obj.K33, obj.SurfaceNewmannData, obj.SurfaceBoundaryEdgeType );
            [ obj.RHS, obj.NonhydrostaticSolver.GlobalStiffMatrix ] = mxAssembleGlobalStiffMatrixWithBottomBCsImposed(...
                obj.NonhydrostaticSolver.GlobalStiffMatrix, obj.RHS, obj.BottomDirichletData,...
                struct(obj.meshUnion.BottomBoundaryEdge), struct(obj.meshUnion.cell), struct(obj.meshUnion), obj.K13, ...
                obj.K23, obj.K33, obj.BottomNewmannData, obj.BottomBoundaryEdgeType );           
            warning('on');
            disp("============For stiff matrix================")
            disp("The maximum difference is:")
            disp(max(max(obj.StiffMatrix - obj.NonhydrostaticSolver.GlobalStiffMatrix)));
            disp("The minimum difference is:")
            disp(min(min(obj.StiffMatrix - obj.NonhydrostaticSolver.GlobalStiffMatrix)));            
            disp("============End stiff matrix================")                   
            obj.SimulatedSolution = obj.NonhydrostaticSolver.GlobalStiffMatrix\obj.RHS(:);
            disp(condest(obj.NonhydrostaticSolver.GlobalStiffMatrix));
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
        
        %> set initial function
        function [fphys2d, fphys] = setInitialField( obj )
            fphys2d = cell( obj.Nmesh, 1 );
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                fphys2d{m} = zeros( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K, obj.Nfield2d );
                fphys{m} = zeros( obj.meshUnion(m).cell.Np, obj.meshUnion(m).K, obj.Nfield );
                fphys2d{m}(:,:,1) = obj.Depth*ones( obj.mesh2d(m).cell.Np, obj.mesh2d(m).K);
            end
        end
        
        function matGetFunction(obj)
            syms x y z nx ny nz;
            obj.Cexact = sin(pi*x)*+sin(pi*y)+sin(pi*z);
%             obj.Cexact = cos(8*pi*x) + cos(8*pi*y) + cos(-8*pi*z);
%             obj.Cexact = sin(-pi/2*z);
            obj.SecondDiffCexact = diff(obj.TempK11*diff(obj.Cexact,x) + obj.TempK13*diff(obj.Cexact, z),x) + ...
                diff(obj.TempK22*diff(obj.Cexact,y) + obj.TempK23*diff(obj.Cexact, z),y) + ...
                diff(obj.TempK13*diff(obj.Cexact,x) + obj.TempK23*diff(obj.Cexact,y) + obj.TempK33*diff(obj.Cexact, z),z);
            
            obj.NewmannCexact = (obj.TempK11*diff(obj.Cexact,x) + obj.TempK13*diff(obj.Cexact, z))*nx + ...
                (obj.TempK22*diff(obj.Cexact,y) + obj.TempK23*diff(obj.Cexact, z))*ny + ...
                (obj.TempK13*diff(obj.Cexact,x) + obj.TempK23*diff(obj.Cexact,y) + obj.TempK33*diff(obj.Cexact, z))*nz;
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

% mesh2d = makeUniformQuadMesh( N, ...
%     [ -obj.ChLength/2, obj.ChLength/2 ], [ -obj.ChWidth/2, obj.ChWidth/2 ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

mesh2d = makeUniformQuadMesh( N, ...
    [ 0, obj.ChLength ], [ 0, obj.ChWidth ], ceil(obj.ChLength/M), ceil(obj.ChWidth/M), bctype);

cell = StdPrismQuad( N, Nz );
zs = zeros(mesh2d.Nv, 1); zb = zs - 1;
% zs = ones(mesh2d.Nv, 1); zb = zs - 1;
mesh3d = NdgExtendMesh3d( cell, mesh2d, zs, zb, Mz );
mesh3d.InnerEdge = NdgSideEdge3d( mesh3d, 1, Mz );
mesh3d.BottomEdge = NdgBottomInnerEdge3d( mesh3d, 1 );
mesh3d.BoundaryEdge = NdgHaloEdge3d( mesh3d, 1, Mz );
mesh3d.BottomBoundaryEdge = NdgBottomHaloEdge3d( mesh3d, 1 );
mesh3d.SurfaceBoundaryEdge = NdgSurfaceHaloEdge3d( mesh3d, 1 );

% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );

end