classdef NdgSideEdge3d < handle
    properties %( SetAccess = protected )
        %> std cell of the three dimensional side edge
        cell
        %> Vandmonde matrix corresponding to the horizontal interpolation point
        V1d
        %> Vandmonde matrix corresponding to the facial interpolation point
        V2d
        %> number of vertical layers
        Nz
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
        %> global interp node index of 1st ele on each edge
        GFToN1
        %> interp node index of 2nd ele on each edge
        FToN2
        %> global interp node index of 2nd ele on each edge
        GFToN2
        %> outward normal vector
        nx, ny, nz
        %> determination of edge Jacabian
        Js
        %> determination of Jacobian in vertical direction
        Jz
        %> length, area or volume of the studied edge
        LAV
    end
    
    methods ( Hidden, Access=protected )
        GetCellSize( obj, mesh )
    end
    
    methods ( Access = public )
        function obj = NdgSideEdge3d( meshUnion3d, meshId, Nz )
            mesh = meshUnion3d( meshId );
            
            obj.Nz = Nz;
            obj.mesh = mesh;
            obj = assembleMassMatrix( obj, mesh.cell.N, mesh.cell.Nz );
            obj = assembleEdgeConnect( obj, mesh );
            obj = assembleNodeProject( obj, mesh );
            obj.Jz = mesh.Jz(obj.GFToN1);
            obj.GetCellSize( mesh );
        end
        
        %> access boundary values at edges
        [ fM, fP ] = matEvaluateSurfValue( obj, fphys );
        %> evaluate strong-form surface term rhs
        [ frhs ] = matEvaluateStrongFromEdgeRHS( obj, fluxM, fluxP, fluxS )
        
        [ value2d ] = VerticalColumnIntegralField( obj, fphys );
    end
    
    methods ( Access = protected )
        obj = assembleEdgeConnect( obj, mesh );
        obj = assembleNodeProject( obj, mesh );
        [ Nfp, M ] = assembleMassMatrix( obj, N, Nz );
        
        [ nx, ny, nz, Js ] = PrismQuadJacobian3d( obj, mesh, f1, e1, fid );
        [ nx, ny, nz, Js ] = PrismTriJacobian3d( obj, mesh, f1, e1, fid );
    end
end