classdef NdgInnerEdge1d < NdgInnerEdge
    %NDGINNEREDGE1D �˴���ʾ�йش����ժҪ
    %   �˴���ʾ��ϸ˵��
    
    properties
    end
    
    methods
        function obj = NdgInnerEdge1d( meshUnion, meshId )
            obj = obj@NdgInnerEdge( meshUnion, meshId );
        end
        
        [ fnode ] = proj_vert2node( obj, fvert );
        
    end
    
    methods ( Access = protected )
        %> connect edge to elements
        obj = assembleMassMatrix( obj );
        obj = assembleEdgeConnect( obj, mesh )
        obj = assembleNodeProject( obj, mesh )
    end
    
end

