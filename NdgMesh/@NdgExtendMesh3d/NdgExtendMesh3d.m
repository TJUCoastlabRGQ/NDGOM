%> @brief three-dimensional vertical extended mesh
classdef NdgExtendMesh3d < handle
    
    properties( Constant )
        type = enumMeshDim.Three
    end
    
    properties ( SetAccess = public )
        %> std cell object
        cell
        %> 2D mesh object
        mesh2d
        %> number of element
        K
        %> vertical layer
        Nz
        %> number of vertices
        Nv
        %> vertex index in each element
        EToV
        %> region types for each element
        EToR
        %> coordinate of vertex
        vx, vy, vz
        %> edge objects
        BottomEdge % surface/bottom faces
        InnerEdge % inner edge
        BottomBoundaryEdge
        BoundaryEdge % halo edge
        SurfaceBoundaryEdge % surface edge
        %> mesh index
        ind
    end
    
    properties
        %> index of adjacent elements and faces
        EToE, EToF
    end
    
    properties
        %> Vertex point to node point relation, data is organized as Node
        %> points number first, followed by global point index
        VertToNode
        %> Weight of each node, and the sum is equal to 1.0
        NodeWeight
    end
    
    % elemental volume infomation
    properties ( SetAccess=protected )
        %> element index in each layer
        EToL
        %> index of adjecent mesh
        EToM
        %> coordinate of interp points
        x, y, z
        %> determination of Jacobian matrix at each interp points
        J, Jz
        %> Jacobian matrix entries
        rx, ry, rz
        sx, sy, sz
        tx, ty, tz
        %> length/area/volume of each element
        LAV
        %> character length of each element
        charLength
        figure_handle
    end
    
    properties ( SetAccess = protected )
        %> central coordinate
        xc, yc, zc
        %> normal vector of each facial point
        %nx, ny, nz
        %> determination of facial integral at each face points
        %Js
    end
    
    methods (Access = public)
        %> construction function
        %> \param[in] cell reference element
        %> \param[in] mesh2d two-dimensional mesh object
        %> \param[in] vzs vertex surface elevation
        %> \param[in] vzb vertex bottom elevation
        %> \param[in] Nz num of vertical layer
        function obj = NdgExtendMesh3d(cell, mesh2d, vzs, vzb, Nz)
            % check vertical elevation
            if( numel(vzs) ~= mesh2d.Nv || numel(vzs) ~= mesh2d.Nv )
                error( 'Input elevation is incorrect.' );
            end
            Kloc = mesh2d.K * Nz;
            Nvert = mesh2d.Nv * (Nz + 1);
            
            obj.cell = cell;
            obj.mesh2d = mesh2d;
            obj.K = Kloc;
            obj.Nv = Nvert;
            obj.Nz = Nz;
            % extend vertex
            obj.vx = repmat(mesh2d.vx, 1, Nz + 1);
            obj.vy = repmat(mesh2d.vy, 1, Nz + 1);
            sigma = linspace(1, 0, Nz + 1);
            obj.vz = vzs * sigma + vzb * (1 - sigma);
            % extend element
            obj.ind = mesh2d.ind;
            obj = ExtendMesh3d( obj, cell, mesh2d, Nz );
            % interp nodes
            obj = GetNodeCoor( obj );
            
            
            % Jacobian coefficients at IP nodes
            [ obj.rx, obj.ry, obj.rz, obj.sx, obj.sy, obj.sz, obj.tx, obj.ty,...
                obj.tz, obj.J, obj.Jz ] = cell.assembleJacobianMatrix( obj.x, obj.y, obj.z);
            obj = GetCellSize( obj );
            
            [ obj.xc ] =  obj.GetMeshAverageValue( obj.x );
            [ obj.yc ] =  obj.GetMeshAverageValue( obj.y );
            [ obj.zc ] =  obj.GetMeshAverageValue( obj.z );
            obj.matVertNodeProjection();
        end% func
        
        function integralValue = GetMeshIntegralValue(obj, nodeVal)
            integralValue = mxGetMeshIntegralValue(nodeVal, obj.cell.wq, obj.J, obj.cell.Vq);
        end
        
        function avergeValue = GetMeshAverageValue(obj, nodeVal)
            integralValue = mxGetMeshIntegralValue(nodeVal, obj.cell.wq, obj.J, obj.cell.Vq);
            avergeValue = integralValue ./ obj.LAV;
        end
        
        function NodeWeightData = GetVertexWeightData( obj, Data)
            NodeWeightData = mxGetVertexWeightedData( obj.VertToNode, obj.NodeWeight, Data );
        end
        
        function nodeQ = proj_vert2node(obj, vertQ)
            % project scalars from mesh verts to nodes
            ele_vQ = vertQ(obj.EToV);
            nodeQ = obj.cell.project_vert2node(ele_vQ);
        end
        %> draw horizontal result
        drawHorizonSlice( obj, field3d )
        %> draw vertical section
        drawVerticalSlice( obj, nodeId2d, field3d )
        
        %> evaluate vertical integral result
        field2d = VerticalColumnIntegralField( obj, field3d )
        %> extend 2d field to 3d
        field3d = Extend2dField( obj, field2d );
        %> evaluate vertical integral from bottom to each node
        fint3d = VerticalIntegralField( obj, field3d )
        
    end
    
    methods ( Access = protected )
        function obj = GetNodeCoor( obj )
            obj.x = obj.proj_vert2node(obj.vx);
            obj.y = obj.proj_vert2node(obj.vy);
            obj.z = obj.proj_vert2node(obj.vz);
        end% func
    end
    
end
