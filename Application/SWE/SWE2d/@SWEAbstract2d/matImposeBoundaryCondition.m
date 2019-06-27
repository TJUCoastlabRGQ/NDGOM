function [ fM, fP ] = matImposeBoundaryCondition( obj, edge, nx, ny, fM, fP, fext )

% impose boundary condition
[ fP ] = mxImposeBoundaryCondition( obj.gra, nx, ny, fP, fext, int8( edge.ftype ) );

fP(:,:,6) = fM(:,:,6);

% hydrostatic reconstruction
[ fM, fP ] = mxHydrostaticReconstruction( obj.hmin, fM, fP );

end

