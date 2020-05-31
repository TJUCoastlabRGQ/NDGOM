function obj = GetCellSize( obj )

% Jacobian determination on each quadrature points
one = ones( obj.cell.Np, obj.K );
% Volume for each cell of the studied mesh
obj.LAV = mxGetMeshIntegralValue( one, obj.cell.wq, obj.J, obj.cell.Vq );

end% func