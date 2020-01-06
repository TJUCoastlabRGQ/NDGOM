function f = EvaluateHorizontalOrthogonalFunc(obj, N, n, r, s)

f = orthogonal_func(obj.Nph, N, n, r, s);

end

function [ fval ] = orthogonal_func(Nph, N, ind, r, s)
ind1 = mod(ind-1, Nph)+1;
j = floor( (ind1 - 1) / (N + 1) );
i = mod( ind1 - 1, N + 1 );

disp('===============');
display(i);
display(j);
disp('===============');

fval = JacobiP( r(:), 0, 0, i ).*JacobiP( s(:), 0, 0, j );

% i = mod( ind - 1, N + 1 );
% j = floor( (ind - 1) / (N + 1) );
% fval = JacobiP( r(:), 0, 0, i ).*JacobiP( s(:), 0, 0, j );
end