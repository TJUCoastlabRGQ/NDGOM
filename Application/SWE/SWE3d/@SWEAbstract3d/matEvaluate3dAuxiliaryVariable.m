function fphys3d = matEvaluate3dAuxiliaryVariable( obj, mesh3d, fphys2d, fphys3d )
%> @brief Function used to calculate the three dimensional auxiliary variable
%> 
%> More detailed description.
%>
%> @param mesh3d The three dimensional mesh object
%> @param fphys3d The three dimensional physical field
%>
%> @retval fphys3d The three dimensional physical field with the auxiliary variable updated

%the vertical unit vector is -1
edge3d = mesh3d.BottomEdge;
[ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
FluxM_1(:, :, 1) = fm(:, :, 1);
FluxM_1(:, :, 2) = fm(:, :, 2);
FluxP_1(:, :, 1) = fp(:, :, 1);
FluxP_1(:, :, 2) = fp(:, :, 2);

%Central flux for the auxiliary variable
fphys3d{1}(:, :, 4:5) =  edge3d.matEvaluateStrongFormEdgeCentralRHS( FluxM_1, FluxP_1 );


%the vertical unit vector is -1
edge3d = mesh3d.BottomBoundaryEdge;
[ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = fm(:, :, 1);
FluxM(:, :, 2) = fm(:, :, 2);
FluxS(:, :, 1) = fp(:, :, 1);   %Zero grad for the primitive variable
FluxS(:, :, 2) = fp(:, :, 2);
fphys3d{1}(:, :, 4:5) = fphys3d{1}(:, :, 4:5) ...
    + edge3d.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );

%the vertical unit vector is 1
edge3d = mesh3d.SurfaceBoundaryEdge;
[ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
FluxM(:, :, 1) = fm(:, :, 1);
FluxM(:, :, 2) = fm(:, :, 2);
FluxS(:, :, 1) = fp(:, :, 1);   %Zero grad for the primitive variable
FluxS(:, :, 2) = fp(:, :, 2);
fphys3d{1}(:, :, 4:5) = fphys3d{1}(:, :, 4:5) ...
    - edge3d.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );

fphys3d{1}(:, :, 6) = mesh3d.Extend2dField( fphys2d{1}(:, :, 4) );
fphys3d{1}(:, :, 7) = mesh3d.Extend2dField( fphys2d{1}(:, :, 1) );

fphys3d{1}(:, :, 4) =  obj.miu./((fphys3d{1}(:, :, 6)).^2) .* (mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,1) ) + fphys3d{1}(:, :, 4) );

fphys3d{1}(:, :, 5) =  obj.miu./((fphys3d{1}(:, :, 6)).^2) .* (mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,2) ) + fphys3d{1}(:, :, 5) );

end