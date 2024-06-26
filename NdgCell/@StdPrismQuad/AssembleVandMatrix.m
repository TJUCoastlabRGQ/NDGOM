function AssembleVandMatrix(obj)
V = zeros(obj.Np, obj.Np);
tempVCV = zeros(obj.Nph, obj.Np);

for n = 1:obj.Np
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r, obj.s );
    fv = obj.EvaluateVerticalOrthogonalFunc( n, obj.t );
    V(:, n) = fh .* fv;
end% for

for n = 1:obj.Np
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r1, obj.s1 );
    fv = obj.EvaluateVerticalOrthogonalFunc( n, zeros(size(obj.r1)) );
    tempVCV(:, n) = fh .* fv;
end

Vh = zeros(obj.Nph, obj.Nph);
for n = 1:obj.Nph
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r1, obj.s1 );
    Vh(:, n) = fh;
end% for

V1d = zeros(obj.Nz+1, obj.Nz+1);
for n = 1:obj.Nz+1
    V1d(:, n) = JacobiP( obj.t1, 0, 0, n-1 );
end% for
% vertical integral vandermonde matrix
[ Vint, VintU ]  = EvaluateVerticalIntegralOrthogonalFunc( obj );

obj.V = V;
obj.V1d = V1d;
obj.Vh = Vh;
obj.Vint = Vint;
obj.VintU = VintU;
obj.VCV = tempVCV/V;
end

% function Vz = EvaluateVerticalIntegralOrthogonalFunc( obj )
% 
% Vz = zeros(obj.Np, obj.Np);
% for n = 1:obj.Np
%     fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r, obj.s );
%     ind = ceil( n / obj.Nph ); % vertical orthogonal polynomial index
%     if ind == 1
%         fv = LegendreNorm1d( ind+1 , obj.t );
%         Vz(:, n) = ( fv - fv(1) ) ./ sqrt( 2 * ind + 1 );
%     else
%         fvm = LegendreNorm1d( ind - 1, obj.t );
%         fvp = LegendreNorm1d( ind + 1, obj.t );
%         Vz(:, n) = ( fvp ./ sqrt( 2*ind+1 ) - fvm ./ sqrt( 2*ind-3 ) ) ./ sqrt( 2*ind-1 );
%     end
%     
%     Vz(:, n) = Vz(:, n) .* fh;
% end% for
% 
% end

function [ Vz, VzU ]  = EvaluateVerticalIntegralOrthogonalFunc( obj )

Vz = zeros(obj.Np, obj.Np);
VzU = zeros(obj.Np, obj.Np);
for n = 1:obj.Np
    fh = obj.EvaluateHorizontalOrthogonalFunc( obj.N, n, obj.r, obj.s );
    ind = ceil( n / obj.Nph ); % vertical orthogonal polynomial index
    if ind == 1
%         fv = LegendreNorm1d( ind+1 , obj.t );
%         Vz(:, n) = ( fv - fv(1) ) ./ sqrt( 2 * ind + 1 );
        fv = LegendreNorm1d( ind , obj.t );
        Vz(:, n) = (obj.t+1) .* fv;
        VzU(:, n) = 2 .* fv - Vz(:,n);
    else
        fvm = LegendreNorm1d( ind - 1, obj.t );
        fvp = LegendreNorm1d( ind + 1, obj.t );
        Vz(:, n) = ( fvp ./ sqrt( 2*ind+1 ) - fvm ./ sqrt( 2*ind-3 ) ) ./ sqrt( 2*ind-1 );
        
        fvm = LegendreNorm1d( ind - 1, ones(size(obj.t)) );
        fvp = LegendreNorm1d( ind + 1, ones(size(obj.t)) );
        VzU(:, n) = ( fvp ./ sqrt( 2*ind+1 ) - fvm ./ sqrt( 2*ind-3 ) ) ./ sqrt( 2*ind-1 ) - ...
            Vz(:, n);
    end
    
    Vz(:, n) = Vz(:, n) .* fh;
    VzU(:, n) = VzU(:, n) .* fh;
end% for

end

function [ f ] = LegendreNorm1d( ind, r )
f = JacobiP(r, 0, 0, ind-1);
end% func