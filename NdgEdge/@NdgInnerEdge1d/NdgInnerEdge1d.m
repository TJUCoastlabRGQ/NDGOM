classdef NdgInnerEdge1d < NdgInnerEdge
    %NDGINNEREDGE1D 此处显示有关此类的摘要
    %   此处显示详细说明
    
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

