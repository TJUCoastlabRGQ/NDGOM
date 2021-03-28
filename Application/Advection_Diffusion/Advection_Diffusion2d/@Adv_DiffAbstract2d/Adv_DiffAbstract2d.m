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
%             [ uNorm ] = fm(:,:,2) .* nx + fm(:,:,3) .* ny;
%             sign_um = sign( uNorm );
%             fluxS = ( fm(:,:,1) .* ( sign_um + 1 ) ...
%                 + fp(:,:,1) .* ( 1 - sign_um  ) ) .* uNorm .* 0.5;
            fluxS = zeros(edge.Nfp, edge.Ne);
            ind = ( edge.ftype == enumBoundaryCondition.Inner );
            [ uNorm ] = fm(:,ind,2) .* nx(:,ind) + fm(:,ind,3) .* ny(:,ind);
            sign_um = sign( uNorm );
            fluxS(:,ind) = ( fm(:,ind,1) .* ( sign_um + 1 ) ...
                + fp(:,ind,1) .* ( 1 - sign_um  ) ) .* uNorm .* 0.5;
            ind = ( edge.ftype == enumBoundaryCondition.Dirichlet );
            [ uNorm ] = fp(:,ind,2) .* nx(:,ind) + fp(:,ind,3) .* ny(:,ind);
            fluxS(:,ind) = fp(:,ind,1) .* uNorm;
            ind = ( edge.ftype == enumBoundaryCondition.Newmann );
            [ uNorm ] = fm(:,ind,2) .* nx(:,ind) + fm(:,ind,3) .* ny(:,ind);
            fluxS(:,ind) = fm(:,ind,1) .* uNorm;
        end
        
        function [ flux ] = matEvaluateSurfFlux( obj, mesh, nx, ny, fm, edge )
            Em = fm(:,:,1) .* fm(:,:,2);
            Gm = fm(:,:,1) .* fm(:,:,3);
            flux = Em .* nx + Gm .* ny;
        end
        
        function [ fm, fp ] = matImposeBoundaryCondition( obj, edge, nx, ny, fm, fp, fext )
            ind = ( edge.ftype == enumBoundaryCondition.Newmann );
            fp(:, ind, 1) = fm(:, ind, 1);
            fp(:, ind, 2) = fm(:, ind, 2);
            fp(:, ind, 3) = fm(:, ind, 3);
            ind = ( edge.ftype == enumBoundaryCondition.Dirichlet );
            fp(:, ind, 1) = fext(:, ind, 1);
            fp(:, ind, 2) = fext(:, ind, 2);
            fp(:, ind, 3) = fext(:, ind, 3);
        end
    end
    
end