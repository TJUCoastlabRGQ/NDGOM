classdef NdgSideEdge3d < handle
    properties ( SetAccess = protected )
        %> type of the side edge
        type
        %> mesh obj
        mesh
        %> num of face nodes
        Nfp
        %> num of edges
        Ne
        %> mass matrix of edge
        M
        %> vertex index on each edge
        FToV
        %> local and adjacent cell index
        FToE
        %> local face index of local and adjacent cell
        FToF
        %> face to mesh index
        FToM
        %> interp node index of 1st ele on each edge
        FToN1
        %> Global interp node index of 1st ele on each edge
        GFToN1
        %> interp node index of 2nd ele on each edge
        FToN2
        %> Global interp node index of 2nd ele on each edge
        GFToN2   
        %> outward normal vector
        nx, ny, nz
        %> determination of edge Jacabian
        Js
    end

    methods ( Access = public )
        function obj = NdgSideEdge3d( meshUnion3d, meshId )
            mesh = meshUnion3d( meshId );
            
            obj.mesh = mesh;
            obj.type = enumStdCell.Quad;
            obj = assembleMassMatrix( obj, mesh.cell.N, mesh.cell.Nz );
            obj = assembleEdgeConnect( obj, mesh );
            obj = assembleNodeProject( obj, mesh );
            obj.assembleGlobalFacialPointIndex;
        end

        %> access boundary values at edges
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );
        %> evaluate strong-form surface term rhs
        [ frhs ] = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS )
    end

    methods ( Access = protected )
        obj = assembleEdgeConnect( obj, mesh );
        obj = assembleNodeProject( obj, mesh );
        [ Nfp, M ] = assembleMassMatrix( obj, N, Nz );

        [ nx, ny, nz, Js ] = PrismQuadJacobian3d( obj, mesh, f1, e1, fid );
        [ nx, ny, nz, Js ] = PrismTriJacobian3d( obj, mesh, f1, e1, fid );
    end
    
    methods( Access = protected )
        function assembleGlobalFacialPointIndex(obj)
            obj.GFToN1 = zeros(size(obj.FToN1));
            obj.GFToN2 = zeros(size(obj.FToN2));
            for i = 1:obj.Ne
               obj.GFToN1(:,i) = ( obj.FToE(1,i) - 1 ) * obj.mesh.cell.Np + ...
                   obj.FToN1(:,i);
               obj.GFToN2(:,i) = ( obj.FToE(2,i) - 1 ) * obj.mesh.cell.Np + ...
                   obj.FToN2(:,i);               
            end
        end
    end    
    
end