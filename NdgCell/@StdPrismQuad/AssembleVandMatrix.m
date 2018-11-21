function AssembleVandMatrix(obj)
V = zeros(obj.Np, obj.Np);

for n = 1:obj.Np
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r, obj.s );
    fv = obj.EvaluateVerticalOrthogonalFunc( n, obj.t );
    V(:, n) = fh .* fv;
end% for

Vh = zeros(obj.Nph, obj.Nph);
for n = 1:obj.Nph
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r1, obj.s1 );
    Vh(:, n) = fh;
end% for

obj.V = V;
obj.Vh = Vh;
end