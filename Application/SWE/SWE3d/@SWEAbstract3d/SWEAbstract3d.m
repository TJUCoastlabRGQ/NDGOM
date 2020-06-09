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
        VerticalEddyViscositySolver
        %>
        HorizontalEddyViscositySolver
    end
    
    properties
        SurfBoundNewmannDate
        BotBoundNewmannDate
        Cf
    end
    
    properties ( SetAccess = public )
        InnerEdgefm3d
        BoundaryEdgefm3d
        InnerEdgefp3d
        BoundaryEdgefp3d
        
        InnerEdgeFluxS3d
        BoundaryEdgeFluxS3d
        InnerEdgeFluxM3d
        BoundaryEdgeFluxM3d
        InnerEdgeFluxP3d
        
        InnerEdgeFluxS2d
        BoundaryEdgeFluxS2d
        InnerEdgeFluxM2d
        BoundaryEdgeFluxM2d
        InnerEdgeFluxP2d
        
        ImplicitRHS3d
    end
    
    properties ( SetAccess = protected )
        %These variable are needed when time stepping
        ExplicitRHS2d = []
        ExplicitRHS3d = []
%         ImplicitRHS3d = []
        % This parameter is used to calculate the horizontal viscosity when smagorinsky model is adopted 
        SmagorinskyConstant = []
        Prantl = []
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
                                                        
        [ fphys3d ] = matEvaluateVerticalVelocity( obj, mesh3d, fphys2d, fphys );
                
        matEvaluateSourceTerm( obj, fphys );
        
        matInitEddyViscositySolver( obj );
        
        matCalculateExplicitRHSTerm( obj, fphys2d, fphys, Stage, RKIndex);
        
        SystemRHS = matAssembleSystemRHS( obj, Tempfphys, SystemRHS, EXa, IMa, dt);
    end
    
    methods ( Sealed, Access = protected )
        %> determine time interval
        [ dt ] = matUpdateTimeInterval( obj, fphys )
    end
    
end