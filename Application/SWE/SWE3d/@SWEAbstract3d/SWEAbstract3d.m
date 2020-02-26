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
        hcrit
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
        %> linear slip parameter
        K
        %> output file object
    end
    
    properties
        %> solver for coriolis source term
        coriolisSolver
        %> solver for unmerical flux
        numfluxSolver
        %> limiter type
        limiterSolver
        %> Solver for the primal continuity equation
        PCESolver2d
        %> Solver for vertical eddy viscosity
        EddyViscositySolver
    end
    
    properties
        WindTaux
        WindTauy
        Cf
    end
    
    methods
        function obj = SWEAbstract3d(  )
            % Doing nothing
        end
        
        initPhysFromOptions( obj, mesh2d, mesh3d );
        AnimationSurfaceLevel( obj );
        
        function matUpdateWetDryState(obj, fphys)
            for m = 1:obj.Nmesh
                wetflag = all( fphys{m}(:,:,1) > obj.hcrit );
                obj.mesh2d(m).status( ~wetflag ) = int8( enumSWERegion.Dry );
                obj.mesh2d(m).status(  wetflag ) = int8( enumSWERegion.Wet );
            end
        end        
    end
    
    % ======================================================================
    methods ( Hidden, Abstract ) % Abstract function, hidden
        %> abstract function to evaluate volume flux term
        [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys );
    end
    % ======================================================================
    
    
    methods( Access = protected )
        
        function matEvaluateTopographySourceTerm( obj, fphys )
            for m = 1:obj.Nmesh
                obj.frhs{m}(:,:,1) = obj.frhs{m}(:,:,1) - obj.gra * fphys{m}(:,:,7) .* fphys{m}(:,:,8);
                obj.frhs{m}(:,:,2) = obj.frhs{m}(:,:,2) - obj.gra * fphys{m}(:,:,7) .* fphys{m}(:,:,9);
            end
        end
        
        function   [ fphys ] = matEvaluatePostFunc(obj, fphys)
              obj.matUpdateWetDryState(fphys);            
        end
                
    end
    
    methods ( Access = protected )
        
        matEvaluateRHS( obj, fphys2d, fphys );
        
        matUpdateExternalField( obj, time, fphys2d, fphys );
        
        matUpdateOutputResult( obj, time, fphys2d, fphys );
        
        matUpdateFinalResult( obj, time, fphys2d, fphys );
        
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
        
        [ fphys3d ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys );
        
        [ TermX, TermY ] = matEvaluateHorizontalPartialDerivativeTerm(obj, mesh3d, fphys);
        
        matEvaluateSourceTerm( obj, fphys );
        
        EddyViscositySolver = matInitEddyViscositySolver( obj );
        
    end
    
    methods ( Sealed, Access = protected )
        %> determine time interval
        [ dt ] = matUpdateTimeInterval( obj, fphys )
    end
    
end

