function [ cellId, Vg ] = accessGaugePointLocation( obj, xg, yg, zg )

Ng = numel( xg );
xb = obj.x([1, end], :);
rd = zeros( Ng, 1 );
cellId = zeros( Ng, 1 );
for n = 1:Ng
    dx = xb - xg(n);
    id = find( dx(1,:).*dx(2,:) <= 0 );
    temprd = ( xg(n) - xb(1, id) )./( xb(2, id) - xb(1, id) )*obj.cell.LAV - 1;
    cellId(n) = id(abs(temprd)<=1);
    rd(n) = temprd(abs(temprd)<=1);
end

Vg = zeros(Ng, obj.cell.Np);
for n = 1:obj.cell.Np
    Vg(:, n) = obj.cell.orthogonal_func(obj.cell.N, n, rd, 0, 0);
end
Vg = Vg/obj.cell.V;

end

