classdef NdgHaloEdge1d < NdgHaloEdge
    %NDGHALOEDGE1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgHaloEdge1d( meshUnion, locMeshId, BCToV )
            obj = obj@NdgHaloEdge( meshUnion, locMeshId, BCToV );
        end
        
         [ fnode ] = proj_vert2node( obj, fvert );
    end
    
    methods( Access = protected )
        obj = assembleMassMatrix( obj );
        obj = assembleEdgeConnect( obj, mesh );
        obj = assembleNodeProject( obj, mesh );
        obj = assembleBoundaryConnection(obj, BCToV);
    end    
    
end

