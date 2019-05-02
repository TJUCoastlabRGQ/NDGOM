classdef SWEAbstract3d < NdgPhysMat
    %SWEABSTRACT3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> cell array field for RHS
        frhs2d
    end
    
    properties ( Abstract )
        %> field to be put in the output file
        outputFieldOrder2d
    end
    
    properties( SetAccess = protected )
        %> cell array for two dimensional external value fields
        fext2d
        %> cell array for three dimensional external value fields
        fext3d
    end
    
    properties ( Abstract, Constant )
        %> wet/dry depth threshold
        hmin
    end
    
%     properties( SetAccess = protected )
%         %> cell array for physical field variable
%         fphys2d
%     end

    properties
        %> cell array for physical field variable
        fphys2d
    end
    properties( Abstract )
        %> number of physical field
        Nfield2d
        %> number of variable field
        Nvar2d
        %> index of variable in physical field
        varFieldIndex2d
    end
    
    properties ( Constant )
        %> gravity acceleration
        gra = 9.81;
    end
    
    properties ( SetAccess = protected )
        %> num of mesh
        %         Nmesh
        %> horizontal mesh
        mesh2d
        %> vertical extended mesh
        mesh3d
        %> viscosity
        miu
        %> linear slip parameter
        K
        %> output file object
    end
    
    properties
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
        %> Solver for the primal continuity equation
        PCESolver2d
    end
    
    properties
        Taux
        Tauy
        Cf
    end
    
    methods
        function obj = SWEAbstract3d(  )
            % Doing nothing
        end
        
        initPhysFromOptions( obj, mesh2d, mesh3d );
        AnimationSurfaceLevel( obj );
    end
    
    % ======================================================================
    methods ( Hidden, Abstract ) % Abstract function, hidden
        %> abstract function to evaluate volume flux term
        [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys );
    end
    % ======================================================================
    
    
    methods( Access = protected )
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            %doing nothing
        end
    end
    
    methods ( Access = protected )
        
        matEvaluateRHS( obj, fphys );
        
        matUpdateExternalField( obj, time, fphys );
        
        matUpdateOutputResult( obj, time, fphys );
        
        matUpdateFinalResult( obj, time, fphys );
        
        matEvaluate2dHorizonMomentum(obj, mesh3d, fphys);
        
        [ VolumeTerm_rhs2d ] = matEvaluate2dHorizonPCEVolumeTerm( obj, mesh2d );
        
        [ InnerSurface_rhs2d ] = matEvaluate2dHorizonPCEInnerSurfaceTerm( obj, InnerEdge);
        
        [ BoundarySurface_rhs2d ] = matEvaluate2dHorizonPCEBoundaryTerm( obj, BoundaryEdge, fext);
        
        [ fphys3d  ] = matEvaluate3dAuxiliaryVariable(  obj, mesh3d, fphys);
        
        [ SideSurface_rhs3d ]  = matEvaluate3dSideSurfaceTerm( obj, InnerEdge, fphys );
        
        [ HorizontalBoundarySurface_rhs3d ] = matEvaluate3dHorizontalBoundaryTerm( obj, BoundaryEdge, fphys, fext );
        
        [ SurfaceBoundary_rhs3d ] = matEvaluate3dSurfaceBoundaryTerm( obj, Edge, fphys );
        
        [ BottomBoundary_rhs3d ] = matEvaluate3dBottomBoundaryTerm( obj, BottomEdge, fphys);
        
        [ AuxialaryVariableFace_rhs3d ] = matEvaluate3dVerticalAuxialaryVariableFaceTerm( obj, mesh3d, fphys );
        
        [ fphys3d ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys );
        
        [ TermX, TermY ] = matEvaluateHorizontalPartialDerivativeTerm(obj, mesh3d, fphys);
        
        matEvaluateSourceTerm( obj, fphys );
        
    end
    
end

