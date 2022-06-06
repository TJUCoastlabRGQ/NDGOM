function  [nx, Js] = matGetElementFaceNormalDirectionVector( obj, ele, adjacentEle,...
    face, BoundaryEdgeFToF, BoundaryEdgenx,...
    BoundaryEdgeFToE, InnerEdgenx, InnerEdgeFToE, BoundaryEdgeJs, InnerEdgeJs )
%> @brief Function to get the unit vector and facial Jacobian at
%> the studied face
%> @details
%> Function to get the unit vector and facial Jacobian at the studied face
%> @param[in] ele The studied local element
%> @param[in] adjacentEle The studied adjacent element
%> @param[in] face The studied local face
%> @param[in] BoundaryEdgeFToF The face to face topological relation of the boundary edge
%> @param[in] BoundaryEdgenx Projection of the unit direction vector in x direction of the boundary edge
%> @param[in] BoundaryEdgeny Projection of the unit direction vector in y direction of the boundary edge
%> @param[in] BoundaryEdgeFToE The face to element topological relation of the boundary edge
%> @param[in] InnerEdgenx Projection of the unit direction vector in x direction of the inner edge
%> @param[in] InnerEdgeny Projection of the unit direction vector in y direction of the inner edge
%> @param[in] InnerEdgeFToE The face to element topological relation of the inner edge
%> @param[in] BoundaryEdgeJs The facial Jacobian of the boundary edge
%> @param[in] InnerEdgeJs The facial Jacobian of the inner edge
%> @param[out] nx Projection of the unit direction vector in x direction of the studied edge
%> @param[out] Js The facial Jacobian of the studied edge
    [nx, Js] = mxGetElementFaceNormalDirectionVector( ele, adjacentEle, face,  BoundaryEdgeFToF, BoundaryEdgenx, ...
        BoundaryEdgeFToE, InnerEdgenx, InnerEdgeFToE, BoundaryEdgeJs, InnerEdgeJs );
end