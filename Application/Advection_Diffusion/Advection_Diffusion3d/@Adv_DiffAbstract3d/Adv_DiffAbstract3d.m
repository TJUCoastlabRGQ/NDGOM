classdef Adv_DiffAbstract3d < Adv_DiffAbstract
    %ADVABSTRACTVARFLOW3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> Number of physical field
        Nfield = 4
        %> Number of variable field
        Nvar = 1
        %> field index of variable field
        varFieldIndex = 1
        
        outputFieldOrder3d = 1
        
        fieldName = {'C'};
    end
    
    properties( Access = protected )
        u0 = 0
        v0 = 0
        w0 = 0
    end
    
    properties
        mesh2d
        %> vertical extended mesh
        mesh3d
    end
    
    properties
        SurfEdgefext
        
        BotEdgefext        
        
        SurfBoundNewmannDate
        
        BotBoundNewmannDate
    end

    methods
        function obj = Adv_DiffAbstract3d()
            obj = obj@Adv_DiffAbstract();
        end
    end
    
    methods ( Abstract, Access = protected )
        f_ext = getExtFunc( obj, mesh, time );
    end
    
    methods ( Access = protected )
                                
        matEvaluateSourceTerm( obj, fphys );
                
    end
    
    methods ( Hidden )
        
        initPhysFromOptions( obj, mesh2d, mesh3d );
        
        function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys )
            E = fphys(:,:,2) .* fphys(:,:,1);
            G = fphys(:,:,3) .* fphys(:,:,1);
            H = fphys(:,:,4) .* fphys(:,:,1);
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, nz, fm, fp, edge )
            [ uNorm ] = fm(:,:,2) .* nx + fm(:,:,3) .* ny + fm(:,:,4) .* nz;
            sign_um = sign( uNorm );
            fluxS = ( fm(:,:,1) .* ( sign_um + 1 ) ...
                + fp(:,:,1) .* ( 1 - sign_um  ) ) .* uNorm .* 0.5;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, nz, fm, edge )
            Em = fm(:,:,1) .* fm(:,:,2);
            Gm = fm(:,:,1) .* fm(:,:,3);
            Hm = fm(:,:,1) .* fm(:,:,4);
            flux = Em .* nx + Gm .* ny + Hm .* nz;
        end
        
        function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, nz, fm, fp, fext )
                ind = ( edge.ftype == enumBoundaryCondition.Clamped );
                fp(:, ind) = fext(:,ind);
                ind = (edge.ftype == enumBoundaryCondition.UpperSurfaceBoundary);
                fp(:, ind) = fext(:,ind);
                ind = (edge.ftype == enumBoundaryCondition.BottomBoundary);
                fp(:, ind) = fext(:,ind);   
        end
    end
    
end

