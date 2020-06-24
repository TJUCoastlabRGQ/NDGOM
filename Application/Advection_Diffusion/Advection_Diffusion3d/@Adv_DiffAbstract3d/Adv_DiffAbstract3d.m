classdef Adv_DiffAbstract3d < NdgPhysMat
    %ADVABSTRACTVARFLOW3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> Number of physical field
        Nfield = 4
        %> Number of variable field
        Nvar = 1
        %> field index of variable field
        varFieldIndex = 1
    end
    
    properties( Access = protected )
        u0 = 0
        v0 = 0
        w0 = 0
    end
    
    properties( Access = protected )
        miu = 0
    end
    
    properties
        outputFieldOrder3d = 1
    end
    
    properties
        SurfBoundNewmannDate
        BotBoundNewmannDate
        BoundaryEdgeNewmannDate
    end
    
    properties ( SetAccess = public )
        InnerEdgefm3d
        BoundaryEdgefm3d
        InnerEdgefp3d
        BoundaryEdgefp3d
        %for shallow water equation in sigma coordinate, the following four
        %part is not used
        SurfaceBoundaryEdgefm3d
        BottomBoundaryEdgefm3d
        SurfaceBoundaryEdgefp3d
        BottomBoundaryEdgefp3d
        ImplicitRHS3d
    end
    
    properties
        %> Solver for vertical eddy viscosity
        VerticalEddyViscositySolver
        %>
        HorizontalEddyViscositySolver
    end    
    
    methods
        function obj = Adv_DiffAbstract3d()
            obj = obj@NdgPhysMat();
        end
    end
    
    methods ( Abstract, Access = protected )
        f_ext = getExtFunc( obj, mesh, time );
    end
    
    methods ( Hidden )
        function initPhysFromOptions( obj, mesh )
            initPhysFromOptions@NdgPhysMat( obj, mesh );
            %here the viscosity solver is to be added
            
            finalTime = obj.getOption('finalTime');
            for m = 1:obj.Nmesh
                obj.fext{m} = obj.getExtFunc( obj.meshUnion(m), finalTime );
            end
        end
        
        function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys )
            E = fphys(:,:,2) .* fphys(:,:,1);
            G = fphys(:,:,3) .* fphys(:,:,1);
            H = fphys(:,:,4) .* fphys(:,:,1);
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, nz, fm, fp )
            [ uNorm ] = fm(:,:,2) .* nx + fm(:,:,3) .* ny + fm(:,:,4) .* nz;
            sign_um = sign( uNorm );
            fluxS = ( fm(:,:,1) .* ( sign_um + 1 ) ...
                + fp(:,:,1) .* ( 1 - sign_um  ) ) .* uNorm .* 0.5;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, edge, nx, ny, nz, fm )
            Em = fm(:,:,1) .* fm(:,:,2);
            Gm = fm(:,:,1) .* fm(:,:,3);
            Hm = fm(:,:,1) .* fm(:,:,4);
            flux = Em .* nx + Gm .* ny + Hm .* nz;
        end
        
        function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, nz, fm, fp, fext )
%             ind = ( edge.ftype == 5 );
%             fp(:, ind) = 0;
        end
    end
    
end

