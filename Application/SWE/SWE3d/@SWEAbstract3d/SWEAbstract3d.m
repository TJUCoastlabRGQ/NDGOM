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
        %> num of mesh
        Nmesh
        %> physical field
        fphys2d, frhs2d, fext2d
        %> physical field
        fphys3d, frhs3d, fext3d
        %> horizontal mesh
        mesh2d
        %> vertical extended mesh
        mesh3d
        %> viscosity
        miu
        %> linear slip parameter
        K
        %> output file object
        outputFile
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
        
        matSolve( obj );
 
        initPhysFromOptions( obj, mesh2d, mesh3d );        
        
    end
    
    
    methods ( Access = protected )
        
        matEvaluateRHS( obj, fphys2d, fphys3d );
        
        matUpdateExternalField( obj, time, fphys2d, fphys3d );
        
        matUpdateOutputResult( obj, time, fphys2d, fphys3d );
        
        [ fphys2d ]  = matEvaluate2dHorizonMomentum(obj, mesh3d, fphys2d, fphys3d);
        
        [ VolumeTerm_rhs2d ] = matEvaluate2dHorizonPCEVolumeTerm( obj, mesh2d, fphys2d );
        
        [ InnerSurface_rhs2d ] = matEvaluate2dHorizonPCEInnerSurfaceTerm( obj, InnerEdge, fphys2d);
        
        [ BoundarySurface_rhs2d ] = matEvaluate2dHorizonPCEBoundaryTerm( obj, BoundaryEdge, fphys2d, fext);
        
        [ fphys3d  ] = matEvaluate3dAuxiliaryVariable(  obj, mesh3d, fphys2d, fphys3d);

        [ SideSurface_rhs3d ]  = matEvaluate3dSideSurfaceTerm( obj, InnerEdge, fphys3d );
        
        [ HorizontalBoundarySurface_rhs3d ] = matEvaluate3dHorizontalBoundaryTerm( obj, BoundaryEdge, fphys3d, fext );
        
        [ SurfaceBoundary_rhs3d ] = matEvaluate3dSurfaceBoundaryTerm( obj, Edge, fphys3d );
        
        [ BottomBoundary_rhs3d ] = matEvaluate3dBottomBoundaryTerm( obj, BottomEdge, fphys3d);
        
        [ AuxialaryVariableFace_rhs3d ] = matEvaluate3dVerticalAuxialaryVariableFaceTerm( obj, mesh3d, fphys3d );
        
        [ fphys3d ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys3d );
        
    end
    
end

