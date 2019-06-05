function obj = assembleMassMatrix( obj )

cell = StdPoint( obj.mesh.cell.N );
obj.M = cell.M;
obj.Nfp = cell.Np;
obj.r = cell.r;

end