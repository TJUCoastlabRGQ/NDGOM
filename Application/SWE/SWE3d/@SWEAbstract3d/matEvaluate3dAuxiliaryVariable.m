function fphys3d = matEvaluate3dAuxiliaryVariable( obj, mesh3d, fphys3d )
%> @brief Function used to calculate the three dimensional auxiliary variable
%> 
%> More detailed description.
%>
%> @param mesh3d The three dimensional mesh object
%> @param fphys3d The three dimensional physical field
%>
%> @retval fphys3d The three dimensional physical field with the auxiliary variable updated

Hmiu = sqrt( obj.miu );

edge3d = mesh3d.BottomEdge;
[ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
FluxM_1(:, :, 1) = fm(:, :, 1); % normal vector is - 1
FluxM_1(:, :, 2) = fm(:, :, 2);
FluxP_1(:, :, 1) = fp(:, :, 1);
FluxP_1(:, :, 2) = fp(:, :, 2);
% FluxS_1 = 0.5 *( FluxM_1 + FluxP_1 );
% FluxS_1 = FluxP_1;
% fphys3d{1}(:, :, 4:5) = edge3d.matEvaluateStrongFormEdgeRHS( FluxM_1, FluxP_1, FluxS_1 );
% fphys3d{1}(:, :, 4:5) = edge3d.matEvaluateStrongFormEdgeAlterRHS( FluxM_1, FluxP_1 );
fphys3d{1}(:, :, 4:5) = edge3d.matEvaluateStrongFormEdgeCentralRHS( FluxM_1, FluxP_1 );

% edge3d = mesh3d.BottomBoundaryEdge;
% [ fm, fp ] = edge3d.matEvaluateSurfValue( fphys3d );
% FluxM(:, :, 1) = fm(:, :, 1);
% FluxM(:, :, 2) = fm(:, :, 2);
% FluxS(:, :, 1) = fp(:, :, 1);
% FluxS(:, :, 2) = fp(:, :, 2);
% fphys3d{1}(:, :, 4:5) = fphys3d{1}(:, :, 4:5) ...
%     + edge3d.matEvaluateStrongFormEdgeRHS( FluxM, FluxS );

fphys3d{1}(:, :, 4) = - Hmiu .* ( fphys3d{1}(:, :, 4) + ...
    mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,1) ) );

fphys3d{1}(:, :, 5) = - Hmiu .* ( fphys3d{1}(:, :, 5) + ...
    mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,2) ) );

% fphys3d{1}(:, :, 4) = - Hmiu .* mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,1) );
% fphys3d{1}(:, :, 5) = - Hmiu .* mesh3d.tz .* ( mesh3d.cell.Dt * fphys3d{1}(:,:,2) );
end