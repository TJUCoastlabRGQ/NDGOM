function [dr, ds, dt] = derivative_orthogonal_func(obj, N1, N2, td, r, s, t)
    td1 = mod( td-1, obj.Nph ) + 1;
    td2 = ceil( td / obj.Nph );

    [ dr, ds ] = quad_derivative_orthogonal_func( N1, td1, r, s);
    [ ft ] = line_orthogonal_func( td2, t );
    % multiply the vertical polynomial
    dr = dr .* ft;
    ds = ds .* ft;

    [ dt ] = line_derivative_orthogonal_func( td2, t);
    [ frs ] = quad_orthogonal_func( N1, td1, r, s );
    % multiply the horizontal polynomial
    dt = dt .* frs;
end

function [ fdr, fds ] = quad_derivative_orthogonal_func( N, ind, r, s)

% transform the index to two indexes.
i = mod( ind - 1, N + 1 );
j = floor( (ind - 1) / (N + 1) );
% calculate the derivative basis function values.
fdr = GradJacobiP(r(:), 0, 0, i).*JacobiP(s(:), 0, 0, j);
fds = JacobiP(r(:), 0, 0, i).*GradJacobiP(s(:), 0, 0, j);
end

function [ f ] = line_orthogonal_func( ind, r )
    f = JacobiP(r, 0, 0, ind-1);
end% func

function [dr] = line_derivative_orthogonal_func( ind, r )
    dr = GradJacobiP(r, 0, 0, ind-1);
end

function [ fval ] = quad_orthogonal_func( N, ind, r, s)

i = mod( ind - 1, N + 1 );
j = floor( (ind - 1) / (N + 1) );
fval = JacobiP( r(:), 0, 0, i ).*JacobiP( s(:), 0, 0, j );

end