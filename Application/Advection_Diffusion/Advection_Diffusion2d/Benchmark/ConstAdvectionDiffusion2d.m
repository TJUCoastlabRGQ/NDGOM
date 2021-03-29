classdef ConstAdvectionDiffusion2d < Adv_DiffAbstract2d
    %CONSTADVECTIONDIFFUSION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    properties
        M
        N
    end
    
    properties
        ExactFunc
        Difft
        SecondDiffx
        SecondDiffy
        Advx
        Advy
        GradInX
        GradInY
    end
    
    methods
        function obj = ConstAdvectionDiffusion2d( N,  M )
            % setup mesh domain
            [ obj.mesh2d  ] = makeChannelMesh( obj, N, M );
            obj.M = M;
            obj.N = N;
            obj.miu = 0.001;
            obj.u0 = 0;
            obj.v0 = 0;
            % allocate boundary field with mesh obj
            obj.matGetExtFunc;
            obj.initPhysFromOptions( obj.mesh2d );
        end
        
        function matGetExtFunc(obj)
            syms x y t;
            obj.ExactFunc = sin(2*pi*t)*sin(2*pi*x)*sin(2*pi*y);
            obj.Difft = diff(obj.ExactFunc, t);
            obj.Advx = obj.u0 * diff(obj.ExactFunc, x);
            obj.Advy = obj.v0 * diff(obj.ExactFunc, y);
            obj.SecondDiffx = diff(obj.miu*diff(obj.ExactFunc, x),x);
            obj.SecondDiffy = diff(obj.miu*diff(obj.ExactFunc, y),y);
            obj.GradInX = obj.miu*diff(obj.ExactFunc,x);
            obj.GradInY = obj.miu*diff(obj.ExactFunc,y);
        end
    end
    
    methods ( Access = protected )
        
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh = obj.meshUnion(m);
                x = mesh.x;
                y = mesh.y;
                t = 0;
                fphys{m}(:,:,1) = eval(obj.ExactFunc);
                fphys{m}(:,:,2) = obj.u0 .* ones(size(fphys{m}(:,:,1)));
                fphys{m}(:,:,3) = obj.v0 .* ones(size(fphys{m}(:,:,1)));
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            x = obj.mesh2d.BoundaryEdge.xb;
            y = obj.mesh2d.BoundaryEdge.yb;
            t = time;
            obj.fext{1}(:,:,1) = eval(obj.ExactFunc);
            obj.fext{1}(:,:,2) = obj.u0 * ones(size(obj.fext{1}(:,:,1)));
            obj.fext{1}(:,:,3) = obj.v0 * ones(size(obj.fext{1}(:,:,1)));
            obj.GradExt = obj.mesh2d.BoundaryEdge.nx .* eval(obj.GradInX) + ...
                obj.mesh2d.BoundaryEdge.ny .* eval(obj.GradInY);
        end
        
        function matEvaluateSourceTerm( obj, time )
            x = obj.meshUnion.x;
            y = obj.meshUnion.y;
            t = time;
            obj.frhs{1} = obj.frhs{1} + ...
                eval(obj.Difft) + ...
                eval(obj.Advx) + ...
                eval(obj.Advy) - ...
                eval(obj.SecondDiffx) - ...
                eval(obj.SecondDiffy);
        end
        
        function f_ext = getExtFunc( obj, mesh, time )
            f_ext = sin(2*pi*time)*sin(2*pi*mesh.x).*sin(pi*mesh.y);
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
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            option('AdvDiffVerticalDiffusionType') = enumVerticalDiffusion.None;
            dx = (obj.mesh2d.cell.r(2) - obj.mesh2d.cell.r(1))/2*(2/obj.M);
            dthu = min( 1/(2*obj.N+1) *  dx/obj.u0, 1/(2*obj.N+1) * dx^2/obj.miu);
            dthv = min( 1/(2*obj.N+1) *  dx/obj.v0, 1/(2*obj.N+1) * dx^2/obj.miu);
%             dthu = min(  1/(2*obj.N+1) * dx^2/obj.miu);
%             dthv = min(  1/(2*obj.N+1) * dx^2/obj.miu);
            option('timeInterval') = 0.02 * min(dthu, dthv);
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('AdvDiffConstantVerticalDiffusionValue') = obj.miu;
            option('AdvDiffHorizontalDiffusionType') = enumHorizontalDiffusion.Constant;
            option('AdvDiffConstantHorizontalDiffusionValue') = obj.miu;
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