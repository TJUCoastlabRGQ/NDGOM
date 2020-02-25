function [ fluxM ] = matEvaluateSurfFlux( obj, edge, nx, ny, nz, fm )

%> $ ( hu^2 + \frac{1}{2}g(H^2 - Z^2) ) * nx + huv * ny + u\omega * nz$
fluxM(:, :, 1) = ( fm(:,:,1).^2./fm(:,:,4)  + 1/2 * obj.gra * ( fm(:,:,4).^2 - fm(:,:,6).^2 ) ) .* nx +...
      ( fm(:,:,1) .* fm(:,:,2) ./ fm(:,:,4) ) .* ny + fm(:,:,1) .* fm(:,:,3)./fm(:,:,4) .* nz; 
%> $huv * nx + ( hv^2+\frac{1}{2}g(H^2 - Z^2) ) * ny + v\omega * nz$
fluxM(:, :, 2) = ( fm(:,:,1) .* fm(:,:,2) ./ fm(:,:,4) ) .* nx +... 
    ( fm(:,:,2).^2./fm(:,:,4) + 1/2 * obj.gra * ( fm(:,:,4).^2 - fm(:,:,6).^2 ) ) .* ny + fm(:,:,2) .* fm(:,:,3)./fm(:,:,4) .* nz; 

% %> $ ( hu^2 + \frac{1}{2}g(H^2) ) * nx + huv * ny + u\omega * nz$
% fluxM(:, :, 1) = ( fm(:,:,1).^2./fm(:,:,4)  + 1/2 * obj.gra * ( fm(:,:,4).^2 ) ) .* nx +...
%       ( fm(:,:,1) .* fm(:,:,2) ./ fm(:,:,4) ) .* ny + fm(:,:,1) .* fm(:,:,3)./fm(:,:,4) .* nz; 
% %> $huv * nx + ( hv^2+\frac{1}{2}g(H^2) ) * ny + v\omega * nz$
% fluxM(:, :, 2) = ( fm(:,:,1) .* fm(:,:,2) ./ fm(:,:,4) ) .* nx +... 
%     ( fm(:,:,2).^2./fm(:,:,4) + 1/2 * obj.gra * ( fm(:,:,4).^2 ) ) .* ny + fm(:,:,2) .* fm(:,:,3)./fm(:,:,4) .* nz; 



end