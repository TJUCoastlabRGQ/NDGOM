classdef SWEAbstract3d < handle
    %SWEABSTRACT3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Abstract, Constant)
        %> wet/dry depth threshold
        hmin
    end
    
    properties (Constant)
        %> gravity acceleration
        gra = 9.81;
    end    
    
    properties ( SetAccess = protected )
        %> gradient of bottom elevation
        zGrad
    end
    
    properties ( SetAccess = private )
        %> solver for coriolis source term
        coriolisSolver
        %> solver for friction source term
        frictionSolver
        %> solver for wind source term
        windSolver
        %> solver for unmerical flux
        numfluxSolver
        %> limiter type
        limiterSolver
    end
    
    methods
        function obj = SWEAbstract3d( N, Nz, M, Mz )
            % Doing nothing
        end
        
        function matSolve( obj )
            matEvaluateRK45( obj );
        end
        
        initPhysFromOptions( obj, mesh2d, mesh3d );        
        
    end
    
    
    methods ( Access = protected )
        
        matEvaluateRHS( obj, fphys2d, fphys3d );
        
        matUpdateExternalField( obj, time, fphys2d, fphys3d );
        
        matUpdateOutputResult( obj, time, fphys2d, fphys3d );
        
        [ fphys3d ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d );
        
        fphys2d = matEvaluate2dHorizonMomentum(obj, mesh3d, fphys2d, fphys3d);
        
        frhs2d_VolumeTerm = matEvaluate2dHorizonPCEVolumeTerm( obj, mesh2d, fphys2d );
        
        frhs2d_InnerSurfaceTerm = matEvaluate2dHorizonPCEInnerSurfaceTerm( obj, InnerEdge, fphys2d);
        
        frhs2d_BoundarySurfaceTerm = matEvaluate2dHorizonPCEBoundaryTerm( obj, BoundaryEdge, fphys2d, fext);
        
        fphys3d = matEvaluate3dAuxiliaryVariable( obj, mesh3d, fphys3d );
        
    end
    
end

