classdef ConstAdvectionDiffusion2d < Adv_DiffAbstract2d
    %CONSTADVECTIONDIFFUSION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    properties
        M
        N
    end
    
    methods
        function obj = ConstAdvectionDiffusion2d( N,  M )
            % setup mesh domain
            [ obj.mesh2d  ] = makeChannelMesh( obj, N, M );
            obj.M = M;
            obj.N = N;
            obj.miu = 0.001;
            obj.u0 = 1;
            obj.v0 = 1;
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d );            
        end
        
        function fext = matGetExtFunc(obj, time)
            fext = obj.getExtFunc(obj.meshUnion(1), time);
        end
    end

    methods ( Access = protected )
        
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1 : obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m}(:,:,1) = getExtFunc(obj, mesh, 0);
                fphys{m}(:,:,2) = obj.u0 .* ones(size(fphys{m}(:,:,1)));
                fphys{m}(:,:,3) = obj.v0 .* ones(size(fphys{m}(:,:,1)));
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            obj.fext{1}(:,:,1) = sin(2*pi*time)*sin(2*pi*obj.mesh2d.BoundaryEdge.xb).*...
                sin(pi*obj.mesh2d.BoundaryEdge.yb);
            obj.fext{1}(:,:,2) = obj.u0 * ones(size(obj.fext{1}(:,:,1)));
            obj.fext{1}(:,:,3) = obj.v0 * ones(size(obj.fext{1}(:,:,1)));        
%         BotBoundNewmannDate
        end
        
        function matEvaluateSourceTerm( obj, time )
            obj.frhs{1} = obj.frhs{1} + ...
                2*pi*cos(2*pi*time)*sin(2*pi*obj.meshUnion.x).*sin(pi*obj.meshUnion.y) + ...
                2*pi*sin(2*pi*time)*cos(2*pi*obj.meshUnion.x).*sin(pi*obj.meshUnion.y) + ...
                pi*sin(2*pi*time)*sin(2*pi*obj.meshUnion.x).*cos(pi*obj.meshUnion.y) + ...
                5*pi^2*obj.miu*sin(2*pi*time)*sin(2*pi*obj.meshUnion.x).*sin(pi*obj.meshUnion.y);
        end
        
        function f_ext = getExtFunc( obj, mesh, time )
            f_ext = sin(2*pi*time)*sin(2*pi*mesh.x).*sin(pi*mesh.y);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 1;
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
            option('timeInterval') = 0.2 * min(dthu, dthv);
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
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped ];

mesh2d = makeUniformTriMesh( N, ...
    [ -1, 1 ], [ -1, 1 ], M, M, bctype);

% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end