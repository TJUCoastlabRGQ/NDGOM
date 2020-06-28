classdef ConstAdvection2d < Adv_DiffAbstract2d
    %CONSTADVECTIONDIFFUSION3D 此处显示有关此类的摘要
    %   此处显示详细说明
    properties
        M
        N
    end
    
    methods
        function obj = ConstAdvection2d( N,  M )
            % setup mesh domain
            [ obj.mesh2d  ] = makeChannelMesh( obj, N, M );
            obj.M = M;
            obj.N = N;
            obj.u0 = 1;
            obj.v0 = 1;
            % allocate boundary field with mesh obj
            obj.initPhysFromOptions( obj.mesh2d );            
        end
        
    end

    methods ( Access = protected )
        
        %> set initial function
        function [ fphys ] = setInitialField( obj )
            fphys = cell( obj.Nmesh, 1 );
            for m = 1:obj.Nmesh
                mesh = obj.meshUnion(m);
                fphys{m}(:,:,1) = getExtFunc(obj, mesh.x, mesh.y, 0);
                fphys{m}(:,:,2) = obj.u0 .* ones(size(fphys{m}(:,:,1)));
                fphys{m}(:,:,3) = obj.v0 .* ones(size(fphys{m}(:,:,1)));
            end
        end
        
        function matUpdateExternalField( obj, time, fphys )
            obj.BoundaryEdgefp{1}(:,:,1) = obj.getExtFunc( obj.mesh2d.BoundaryEdge.xb, obj.mesh2d.BoundaryEdge.yb, time);
            obj.BoundaryEdgefp{1}(:,:,2) = obj.u0 .* ones(size(obj.BoundaryEdgefp{1}(:,:,1)));
            obj.BoundaryEdgefp{1}(:,:,3) = obj.v0 .* ones(size(obj.BoundaryEdgefp{1}(:,:,1)));
        end
        
        function matEvaluateSourceTerm( obj, time )
            obj.frhs{1} = obj.frhs{1} + ...
                2*pi*cos(2*pi*time)*sin(2*pi*obj.meshUnion.x).*sin(pi*obj.meshUnion.y) + ...
                2*pi*sin(2*pi*time)*cos(2*pi*obj.meshUnion.x).*sin(pi*obj.meshUnion.y) + ...
                pi*sin(2*pi*time)*sin(2*pi*obj.meshUnion.x).*cos(pi*obj.meshUnion.y);
        end
        
        function f_ext = getExtFunc( obj, x, y, time )
            f_ext = sin(2*pi*time)*sin(2*pi*x).*sin(pi*y);
        end
        
        function [ option ] = setOption( obj, option )
            ftime = 1;
            outputIntervalNum = 500;
            option('startTime') = 0.0;
            option('finalTime') = ftime;
            option('outputIntervalType') = enumOutputInterval.DeltaTime;
            option('outputTimeInterval') = ftime/outputIntervalNum;
            option('outputCaseName') = mfilename;
            option('outputNcfileNum') = 1;
            option('temporalDiscreteType') = enumTemporalDiscrete.IMEXRK343;
            dx = (obj.mesh2d.cell.r(2) - obj.mesh2d.cell.r(1))/2*(2/obj.M);
            dthu = 1/(2*obj.N+1) *  dx/obj.u0;
            dthv = 1/(2*obj.N+1) *  dx/obj.v0; 
            option('timeInterval') = 0.6 *  min(dthu, dthv);
            option('equationType') = enumDiscreteEquation.Strong;
            option('integralType') = enumDiscreteIntegral.QuadratureFree;
            option('outputType') = enumOutputFile.NetCDF;
            option('AdvDiffHorizontalDiffusionType') = enumHorizontalDiffusion.None;
        end
        
    end    
end

function [ mesh2d ] = makeChannelMesh( obj, N, M )

bctype = [ ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped, ...
    enumBoundaryCondition.Clamped ];

mesh2d = makeUniformQuadMesh( N, ...
    [ -1, 1 ], [ -1, 1 ], M, M, bctype);

% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'West-East' );
% [ mesh2d, mesh3d ] = ImposePeriodicBoundaryCondition3d(  mesh2d, mesh3d, 'South-North' );
end