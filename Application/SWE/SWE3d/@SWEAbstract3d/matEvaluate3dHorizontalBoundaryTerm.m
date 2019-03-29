function frhs3d = matEvaluate3dHorizontalBoundaryTerm( obj, edge, fphys3d, fext3d )
%EVALUATE3D_HORIZONTALBOUNDARYKERNEL Summary of this function goes here
%   Detailed explanation goes here

[ fm, fp ] = edge.matEvaluateSurfValue( fphys3d );

% apply clamped boundary condition
ind = ( edge.ftype == enumBoundaryCondition.Clamped );
fp(:, ind, 1) = fext3d(:, ind, 1);
fp(:, ind, 2) = fext3d(:, ind, 2);

ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
Hun =  fm( :, ind, 1 ) .* edge.nx(:, ind) + fm( :, ind, 2).* edge.ny(:, ind);
Hvn = -fm( :, ind, 1 ) .* edge.ny(:, ind) + fm( :, ind, 2).* edge.nx(:, ind);

fp(:, ind, 1) = - Hun .* edge.nx(:, ind) - Hvn .* edge.ny(:, ind);
fp(:, ind, 2) = - Hun .* edge.ny(:, ind) + Hvn .* edge.nx(:, ind);

%> $(hu^2+\frac{1}{2}g(H^2 - h^2) )* nx + huv * ny + u\omega * nz$
FluxM(:, :, 1) = ( fm(:,:,1).^2./fm(:,:,6)  + 1/2 * obj.gra * ( fm(:,:,6).^2 - obj.H0^2 ) ) .* edge.nx +...
      ( fm(:,:,1) .* fm(:,:,2) ./ fm(:,:,6) ) .* edge.ny;
%> $(hu^2+\frac{1}{2}g(H^2 - h^2) )* nx + huv * ny + u\omega * nz$
FluxP(:, :, 1) = ( fp(:,:,1).^2./fp(:,:,6) + 1/2 * obj.gra * ( fp(:,:,6).^2 - obj.H0^2 )) .* edge.nx +...
      ( fp(:,:,1) .* fp(:,:,2) ./ fp(:,:,6) ) .* edge.ny;
%> $huv * nx + ( hv^2+\frac{1}{2}g(H^2 - h^2) ) * ny + v\omega * nz$
FluxM(:, :, 2) = ( fm(:,:,1) .* fm(:,:,2) ./ fm(:,:,6) ) .* edge.nx +... 
    ( fm(:,:,2).^2./fm(:,:,6) + 1/2 * obj.gra * ( fm(:,:,6).^2 - obj.H0^2 ) ) .* edge.ny; 
%> $huv * nx + ( hv^2+\frac{1}{2}g(H^2 - h^2) ) * ny + v\omega * nz$
FluxP(:, :, 2) = ( fp(:,:,1) .* fp(:,:,2) ./ fp(:,:,6) ) .* edge.nx +... 
    ( fp(:,:,2).^2./fp(:,:,6) + 1/2 * obj.gra * ( fp(:,:,6).^2 - obj.H0^2 ) ) .* edge.ny; 

%> $\lambda = abs( max( u*nx + v*ny + \frac{\omega}{H}*nz, ( u+\sqrt(gH) )*nx + ( v+\sqrt(gH) )*ny + \frac{\omega}{H}*nz,( u-\sqrt(gH) )*nx + ( v-\sqrt(gH) )*ny + \frac{\omega}{H}*nz )  )$
lambda = abs( max( vertcat (max ( max ( fm(:, :, 1) ./ fm(:, :, 6) .* edge.nx +  fm(:, :, 2) ./ fm(:, :, 6) .* edge.ny ,...
    - fp(:, :, 1) ./ fp(:, :, 6) .* edge.nx -  fp(:, :, 2) ./ fp(:, :, 6) .* edge.ny )), ...
    max( max (( fm(:, :, 1) ./ fm(:, :, 6) + sqrt(obj.gra* fm(:, :, 6))) .* edge.nx...
      + ( fm(:, :, 2) ./ fm(:, :, 6) + sqrt(obj.gra* fm(:, :, 6))) .* edge.ny ,...
   ( -fp(:, :, 1) ./ fp(:, :, 6) - sqrt(obj.gra* fp(:, :, 6))) .* edge.nx...
      + ( -fp(:, :, 2) ./ fp(:, :, 6) - sqrt(obj.gra* fp(:, :, 6))) .* edge.ny )),...
  max( max (( fm(:, :, 1) ./ fm(:, :, 6) - sqrt(obj.gra* fm(:, :, 6))) .* edge.nx...
      + ( fm(:, :, 2) ./ fm(:, :, 6) - sqrt(obj.gra* fm(:, :, 6))) .* edge.ny,...
   ( -fp(:, :, 1) ./ fp(:, :, 6) +  sqrt(obj.gra* fp(:, :, 6))) .* edge.nx...
      + ( -fp(:, :, 2) ./ fp(:, :, 6) + sqrt(obj.gra* fp(:, :, 6))) .* edge.ny)))));
  

%> $\mathbf n\cdot\mathbf {F^*} = \frac{\mathbf{F^{(+)}}+\mathbf{F^{(-)}}}{2} - \frac{\lambda}{2}(\bold U^+ - \bold U^-)$
FluxS(:, :, 1) = 0.5 .* ( FluxM(:, :, 1) + FluxP(:, :, 1) - ...
    bsxfun( @times, lambda , ( fp( :, :, 1 ) - fm( :, :, 1 ) ) ) );
FluxS(:, :, 2) = 0.5 .* ( FluxM(:, :, 2) + FluxP(:, :, 2) - ...
    bsxfun( @times, lambda , ( fp( :, :, 2 ) - fm( :, :, 2 ) ) ) );

frhs3d = edge.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );

end

