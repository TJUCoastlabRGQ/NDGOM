function obj = assembleMassMatrix( obj, N, Nz )

Nph = N + 1;
Npz = Nz + 1;
Nfp = Nph * Npz;

[ r, s ] = assembleInpterpNode( N, Nz );

V = zeros( Nfp, Nfp );
for i = 1 : Nfp
    V(:, i) = orthogonalFunc( i, N, r, s );
end

invV = inv(V);
M = (invV')*invV;

cell = StdLine(N);
obj.V1d = cell.V;
obj.V2d = V;
obj.Nfp = Nfp;
obj.M = M;
end