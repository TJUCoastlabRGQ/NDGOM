classdef EllipticMixedParticalDerivativeTest2d < Adv_DiffAbstract2d
    %ELLIPTICMIXEDPARTICALDERIVATIVETEST2D 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties
        StiffMatrix
        PtProducedStiffMatrix
    end
    
    properties
        ExactFunc
        RHS
        DirichletData
        NewmannData
        ExactSolution
        SimulatedSolution
        MixedSecondDiffTerm
    end
    
    methods
        
        function obj = EllipticMixedParticalDerivativeTest2d(N, M)
            % setup mesh domain
            [ obj.mesh2d  ] = makeChannelMesh( obj, N, M );
            obj.initPhysFromOptions( obj.mesh2d );
            obj.HorizontalEddyViscositySolver = MixedHorzDiffSolver(obj);
            obj.matGetExtFunc;
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            obj.RHS = eval(obj.MixedSecondDiffTerm);
            obj.ExactSolution = eval(obj.ExactFunc);
            obj.AssembleGlobalStiffMatrix;
            obj.matGetStiffMatrixInPointForm;
        end
        
        function EllipticProblemSolve(obj)
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            obj.SimulatedSolution = obj.PtProducedStiffMatrix\obj.RHS(:);
        end
        
        function matGetExtFunc(obj)
            syms x y;
%             obj.ExactFunc = sin(2*pi*x)*sin(pi/2*y);
%             obj.MixedSecondDiffTerm = diff(diff(obj.ExactFunc, y),x);
%             DiffFunc = diff(obj.ExactFunc,y);
%             x = obj.meshUnion.BoundaryEdge.xb;
%             y = obj.meshUnion.BoundaryEdge.yb;
%             obj.NewmannData = obj.meshUnion.BoundaryEdge.nx .* obj.meshUnion.BoundaryEdge.nx .* eval(DiffFunc);
%             obj.DirichletData = eval(obj.ExactFunc);
            obj.ExactFunc = sin(2*pi*x)*sin(2*pi*y);
            obj.MixedSecondDiffTerm = diff(diff(obj.ExactFunc, x),x) + diff(diff(obj.ExactFunc, y),y);
            DiffFunc = diff(obj.ExactFunc,y);
            x = obj.meshUnion.BoundaryEdge.xb;
            y = obj.meshUnion.BoundaryEdge.yb;
            obj.NewmannData = obj.meshUnion.BoundaryEdge.nx .* obj.meshUnion.BoundaryEdge.nx .* eval(DiffFunc) + obj.meshUnion.BoundaryEdge.ny .* obj.meshUnion.BoundaryEdge.ny .* eval(DiffFunc);
            obj.DirichletData = eval(obj.ExactFunc);
        end
        
        AssembleGlobalStiffMatrix(obj);
    end
    
    methods ( Access = protected )
        
        function matGetStiffMatrixInPointForm( obj )
            K = obj.meshUnion.K;
            Np = obj.meshUnion.cell.Np;
            obj.PtProducedStiffMatrix = zeros(K*Np);
            fphys = cell(1);
            obj.fext = cell(1);
            obj.GradExt = zeros(obj.meshUnion.BoundaryEdge.Nfp , obj.meshUnion.BoundaryEdge.Ne);
            for i = 1:obj.meshUnion.K*obj.meshUnion.cell.Np
                fphys{1} = zeros(obj.meshUnion.cell.Np , obj.meshUnion.K);
                fphys{1}(i) = 1;
                obj.PtProducedStiffMatrix((i-1)*K*Np+1:i*K*Np) = obj.HorizontalEddyViscositySolver.matEvaluateStiffMatrixInPointForm( obj, fphys);
            end
        end
        
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m}(:,:,1) = zeros(size(mesh.x));
            end
        end
        
        
        function f_ext = getExtFunc( obj, mesh, time )
            f_ext = sin(2*pi*mesh.x).*sin(pi/2*mesh.y);
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
            option('AdvDiffHorizontalDiffusionType') = enumHorizontalDiffusion.Constant;
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
    enumBoundaryCondition.Dirichlet, ...
    enumBoundaryCondition.Dirichlet, ...
    enumBoundaryCondition.Dirichlet, ...
    enumBoundaryCondition.Dirichlet ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -1, 1 ], [ -1, 1 ], M, M, bctype);

% [ mesh2d ] = ImposePeriodicBoundaryCondition2d(  mesh2d, 'West-East' );
% [ mesh2d ] = ImposePeriodicBoundaryCondition2d(  mesh2d, 'South-North' );
end
