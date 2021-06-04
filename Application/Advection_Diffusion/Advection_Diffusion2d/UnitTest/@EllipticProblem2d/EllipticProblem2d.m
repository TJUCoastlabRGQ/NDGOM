classdef EllipticProblem2d < Adv_DiffAbstract2d
    %ELLIPTICPROBLEM2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        StiffMatrix
    end
    
    properties
        D11
        D12
        D21
        D22
    end
    
    
    properties
        ExactFunc
        RHS
        DirichletData
        NewmannData
        ExactSolution
        SimulatedSolution
        SecondDiffTerm
    end
    
    methods
        
        function obj = EllipticProblem2d(N, M, DiffCoe)
            % setup mesh domain
            [ obj.mesh2d  ] = makeChannelMesh( obj, N, M );
            obj.initPhysFromOptions( obj.mesh2d );
            obj.D11 = DiffCoe(1);
            obj.D12 = DiffCoe(2);
            obj.D21 = DiffCoe(3);
            obj.D22 = DiffCoe(4);
%             obj.HorizontalEddyViscositySolver = MixedHorzDiffInXSolver(obj);
            obj.matGetExtFunc;
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            obj.RHS = eval(obj.SecondDiffTerm);
%             obj.RHS = zeros(size(obj.meshUnion.x));
            obj.ExactSolution = eval(obj.ExactFunc);
            obj.AssembleGlobalStiffMatrix;
        end
        
        function EllipticProblemSolve(obj)
            disp(condest(obj.StiffMatrix));
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            obj.SimulatedSolution = obj.StiffMatrix\obj.RHS(:);
        end
        
        function matGetExtFunc(obj)
            syms x y;
%             obj.ExactFunc = sin(pi*x)*sin(pi/2*y);
            obj.ExactFunc = sin(-pi/2*x) + sin(-pi/2*y);
            obj.SecondDiffTerm = diff(obj.D11 * diff(obj.ExactFunc, x),x) + diff(obj.D12 * diff(obj.ExactFunc, y),x) + ...
                diff(obj.D21 * diff(obj.ExactFunc, x),y) + diff(obj.D22 * diff(obj.ExactFunc, y),y);
            DiffFuncX = diff(obj.ExactFunc, x);
            DiffFuncY = diff(obj.ExactFunc, y);
            x = obj.meshUnion.BoundaryEdge.xb;
            y = obj.meshUnion.BoundaryEdge.yb;
            obj.NewmannData = obj.meshUnion.BoundaryEdge.nx .* eval(obj.D11 * DiffFuncX + obj.D12 * DiffFuncY ) + ...
                obj.meshUnion.BoundaryEdge.ny .* eval(obj.D21 * DiffFuncX + obj.D22 * DiffFuncY );
%             obj.NewmannData = obj.meshUnion.BoundaryEdge.nx .* eval(obj.D11 * DiffFuncX ) + ...
%                 obj.meshUnion.BoundaryEdge.ny .* eval( obj.D22 * DiffFuncY );            
            obj.DirichletData = eval(obj.ExactFunc);
        end
        
        AssembleGlobalStiffMatrix(obj);
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m}(:,:,1) = zeros(size(mesh.x));
            end
        end
        
        
        function f_ext = getExtFunc( obj, mesh, time )
            f_ext = sin(-pi/2*mesh.x) + sin(-pi/2*mesh.y);
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
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('AdvDiffHorizontalDiffusionType') = enumHorizontalDiffusion.None;
            option('AdvDiffConstantHorizontalDiffusionValue') = 1;
        end
    end
    
    methods(Hidden)
        function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, fm, fp, fext )
            ind = ( edge.ftype == enumBoundaryCondition.Newmann );
            fp(:, ind, 1) = fm(:, ind, 1);
            ind = ( edge.ftype == enumBoundaryCondition.Dirichlet );
            fp(:, ind, 1) = 0;
        end
    end
    
end

function [ mesh2d ] = makeChannelMesh( obj, N, M )

bctype = [ ...
    enumBoundaryCondition.Newmann, ...
    enumBoundaryCondition.Newmann, ...
    enumBoundaryCondition.Newmann, ...
    enumBoundaryCondition.Newmann ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -1.8, 1.8 ], [ -1.8, 1.8 ], M, M, bctype);

% [ mesh2d ] = ImposePeriodicBoundaryCondition2d(  mesh2d, 'West-East' );
% [ mesh2d ] = ImposePeriodicBoundaryCondition2d(  mesh2d, 'South-North' );
end

