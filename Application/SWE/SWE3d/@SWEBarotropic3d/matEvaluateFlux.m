function [ E, G, H ] = matEvaluateFlux( obj, mesh, fphys3d )

%> $Hu^2 + \frac{1}{2}g(H^2-h^2)$ at present only constant bathmetry considered
E(:,:,1) = fphys3d(:, :, 1).^2./fphys3d(:, :, 4) + 1/2 * obj.gra * ( fphys3d(:, :, 4).^2 ); 
% E(:,:,1) = fphys3d(:, :, 1).^2./fphys3d(:, :, 4) + 1/2 * obj.gra * ( fphys3d(:, :, 4).^2 - obj.H0^2 ); 
%> $Huv$
G(:,:,1) = fphys3d(:, :, 1) .* fphys3d(:, :, 2)./fphys3d(:, :, 4);
%> $\omega u$
H(:,:,1) = fphys3d(:, :, 1) .* fphys3d(:, :, 3)./fphys3d(:, :, 4); 

%> $Huv$
E(:,:,2) = fphys3d(:, :, 1) .* fphys3d(:, :, 2)./fphys3d(:, :, 4); 
%> $Hv^2 + \frac{1}{2}g(H^2-h^2)$
G(:,:,2) = fphys3d(:, :, 2).^2./fphys3d(:, :, 4) + 1/2 * obj.gra * ( fphys3d(:, :, 4).^2 );
% G(:,:,2) = fphys3d(:, :, 2).^2./fphys3d(:, :, 4) + 1/2 * obj.gra * ( fphys3d(:, :, 4).^2 - obj.H0^2 );
%> $\omega v$
H(:,:,2) = fphys3d(:, :, 2) .* fphys3d(:, :, 3)./fphys3d(:, :, 4);  
end