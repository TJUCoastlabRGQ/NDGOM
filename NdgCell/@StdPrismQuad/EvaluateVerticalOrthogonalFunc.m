function f = EvaluateVerticalOrthogonalFunc( obj, n, t )

td2 = ceil( n / obj.Nph );
[ f ] = line_orthogonal_func( td2, t );

end

function [ f ] = line_orthogonal_func( ind, r )
f = JacobiP(r, 0, 0, ind-1);
end% func
