function GetCellSize( obj, mesh )
% Jacobian determination on each quadrature points
one = ones( obj.cell.Np, obj.Ne );
% Volume for each cell of the studied mesh
obj.LAV = mxGetMeshIntegralValue( one, obj.cell.wq, obj.Js, obj.cell.Vq );
end