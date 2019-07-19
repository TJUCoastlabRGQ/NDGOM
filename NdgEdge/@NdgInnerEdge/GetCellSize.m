function obj = GetCellSize( obj )

% Jacobian determination on each quadrature points
one = ones( obj.cell.Np, obj.Ne );
LAV = mxGetMeshIntegralValue( one, obj.cell.wq, obj.Js, obj.cell.Vq );
obj.LAV = LAV;

end% func