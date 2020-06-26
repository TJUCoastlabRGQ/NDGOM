classdef Adv_DiffAbstract2d< Adv_DiffAbstract
    %ADVABSTRACTVARFLOW3D Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %> Number of physical field
        Nfield = 3
        %> Number of variable field
        Nvar = 1
        %> field index of variable field
        varFieldIndex = 1
        
        fieldName = {'C'};
    end
    
    properties( Access = protected )
        
        u0 = 0
        
        v0 = 0
    end
    
    properties( SetAccess = protected )
        %> cell array for three dimensional external value fields
        fext2d
    end
    
    properties
        
        mesh2d
        
        outputFieldOrder2d = 1
        
        outputFile2d
    end
    
    methods
        function obj = Adv_DiffAbstract2d()
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
        
        
        
        function [ E, G ] = matEvaluateFlux( obj, mesh, fphys )
            E = fphys(:,:,2) .* fphys(:,:,1);
            G = fphys(:,:,3) .* fphys(:,:,1);
        end
        
        function [ fluxS ] = matEvaluateSurfNumFlux( obj, mesh, nx, ny, fm, fp, edge )
            [ uNorm ] = fm(:,:,2) .* nx + fm(:,:,3) .* ny;
            sign_um = sign( uNorm );
            fluxS = ( fm(:,:,1) .* ( sign_um + 1 ) ...
                + fp(:,:,1) .* ( 1 - sign_um  ) ) .* uNorm .* 0.5;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, edge, nx, ny, fm )
            Em = fm(:,:,1) .* fm(:,:,2);
            Gm = fm(:,:,1) .* fm(:,:,3);
            flux = Em .* nx + Gm .* ny;
        end
        
        function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, nz, fm, fp, fext )
            %             ind = ( edge.ftype == 5 );
            %             fp(:, ind) = 0;
        end
    end
    
end