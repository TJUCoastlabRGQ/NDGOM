classdef Adv_DiffAbstract < NdgPhysMat
    %ADVABSTRACTVARFLOW3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties( Abstract )
        %> Number of physical field
        Nfield
        %> Number of variable field
        Nvar
        %> field index of variable field
        varFieldIndex
    end
    
    properties( Access = protected )
        miu = 0

    end
    
    properties
        outputFile
    end
    
    properties ( SetAccess = public )
        
        ExplicitRHS
        ImplicitRHS
        %The exact value at the final time point
        ExactValue        
    end
    
    properties
        %>
        HorizontalEddyViscositySolver
        
        %> Solver for vertical eddy viscosity
        VerticalEddyViscositySolver
    end
    
    methods
        function obj = Adv_DiffAbstract()
            obj = obj@NdgPhysMat();
        end
    end
    
    methods ( Abstract, Access = protected )
        f_ext = getExtFunc( obj, mesh, time );
    end
    
    methods ( Access = protected )
        
        matUpdateOutputResult( obj, time, fphys );
        
        matUpdateFinalResult( obj, time, fphys );
                
    end
    
    methods ( Hidden, Abstract )
        initPhysFromOptions( obj, mesh );
        
        
        [ E, G ] = matEvaluateFlux( obj, mesh, fphys );
        
        
        [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, fm, fp, edge );
        
        
        [ flux ] = matEvaluateSurfFlux( obj, edge, nx, fm );
        
        
        [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx,  fm, fp, fext );
        
    end
    
end