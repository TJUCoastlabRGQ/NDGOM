function [Nq, rq, sq, tq, wq] = quad_coor_func( obj, N, N2 )
    [ Nq1, rq1, sq1, wq1 ] = quad_quadrature_coor_func( N );
    [ Nq2, tq2, wq2 ] = line_quadrature_coor_func( N2 );

    Nq = Nq1 * Nq2;
    r = repmat(rq1, 1, Nq2);
    s = repmat(sq1, 1, Nq2);
    wq1 = repmat(wq1, 1, Nq2);
    t = repmat(tq2', Nq1, 1 );
    wq2 = repmat(wq2', Nq1, 1 );
    wq = wq1 .* wq2;

    rq = r(:); sq = s(:); tq = t(:); wq = wq(:);
end

function [ Nq,rq,wq ] = line_quadrature_coor_func( N )
    Nq = 6*N;
    [ rq, wq ] = zwglj( Nq );
end

function [ Nq, rq, sq, wq ] = quad_quadrature_coor_func( N )
% np = N+1;
np = 6 * N;
% [zq, w] = zwgl(np); % the 1D LGL quadrature points and their weights
[ zq, w ] = zwglj(np);
% loop along the r-axis first
rq = zq*ones(1, np);
sq = ones(np, 1)*zq';
rq = rq(:); 
sq = sq(:); 
wq = bsxfun(@times, w, w'); 
wq = wq(:);
tq = zeros(size(rq));
Nq = numel(rq);
end% func