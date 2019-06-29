function matCalculateLDGPenaltyParameter( obj, mesh)
%> @brief Function to calculate the penalty parameter C11 for the LDG flux calculation
%> @details
%> Function to  calculate the penalty parameter C11 for the LDG flux calculation.
%> This function is programmed according to Hesthaven and Warburton 2007.
%> @param[in] mesh The mesh object
InnerEdge = mesh.InnerEdge;
BoundaryEdge = mesh.BoundaryEdge;
N = mesh.cell.N;
% obj.IETau = max ( InnerEdge.Js ./ mesh.J(InnerEdge.FToN1 +...
%     ( InnerEdge.FToE(1,:) - 1 ) * mesh.cell.Np ),...
%    InnerEdge.Js ./ mesh.J(InnerEdge.FToN2 + ...
%    ( InnerEdge.FToE(2,:) - 1 ) * mesh.cell.Np ) ) * (N + 1) ^2 * 100;
% obj.BETau =  BoundaryEdge.Js ./ mesh.J(BoundaryEdge.FToN1 +...
%     ( BoundaryEdge.FToE(1,:) - 1 ) * mesh.cell.Np ) * (N + 1) ^2 * 100;
obj.IETau = max ( InnerEdge.Js ./ mesh.J( bsxfun(@plus, InnerEdge.FToN1, ...
    ( InnerEdge.FToE(1,:) - 1 ) * mesh.cell.Np )),...
   InnerEdge.Js ./ mesh.J( bsxfun(@plus, InnerEdge.FToN2, ...
    ( InnerEdge.FToE(2,:) - 1 ) * mesh.cell.Np ) ) ) * (N + 1) ^2 * 100;
obj.BETau =  BoundaryEdge.Js ./ mesh.J( bsxfun ( @plus , BoundaryEdge.FToN1, ...
    ( BoundaryEdge.FToE(1,:) - 1 ) * mesh.cell.Np )) * (N + 1) ^2 * 100;
end