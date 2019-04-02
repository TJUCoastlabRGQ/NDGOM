function fphys = matEvaluateCorrectedMomentum( obj, mesh3d, fphys, fphys2d )
%> @brief Function used to correct the three dimensional horizontal momentum
%>
%> More detailed description.
%>
%> @param mesh3d The three dimensional mesh object
%> @param fphys The three dimensional physical field
%> @param fphys2d The two dimensional physical field
%>
%> @retval fphys The three dimensional physical field with horizontal momentum updated
% tempHorMomentum3d = cell(obj.Nmesh);
% data = cell(1);
% for m = 1:obj.Nmesh
%     tempHorMomentum3d{m}(:,:,1) = mesh3d.VerticalIntegralField( fphys{m}(:,:,1) );
%     tempHorMomentum3d{m}(:,:,2) = mesh3d.VerticalIntegralField( fphys{m}(:,:,2) );
%     data(1) = tempHorMomentum3d{m};
%     
%     edge3d = mesh3d(m).SurfaceBoundaryEdge;
%     [ fm, ~ ] = edge3d.matEvaluateSurfValue( data );
%     fphys{m}(:,:,1) = bsxfun( @plus, fphys{m}(:,:,1), fphys2d{m}(:,:,2) - fm(:,:,1) );
%     fphys{m}(:,:,2) = bsxfun( @plus, fphys{m}(:,:,2), fphys2d{m}(:,:,3) - fm(:,:,2) );
% end
    data = cell(1);
    tempHorMomentum3d(:,:,1) = mesh3d.VerticalIntegralField( fphys(:,:,1) );
    tempHorMomentum3d(:,:,2) = mesh3d.VerticalIntegralField( fphys(:,:,2) );
    data{1} = tempHorMomentum3d;
    
    edge3d = mesh3d.SurfaceBoundaryEdge;
    [ fm, ~ ] = edge3d.matEvaluateSurfValue( data );
    fphys(:,:,1) = bsxfun( @plus, fphys(:,:,1), mesh3d.Extend2dField( fphys2d(:,:,2) - fm(:,:,1) ) );
    fphys(:,:,2) = bsxfun( @plus, fphys(:,:,2), mesh3d.Extend2dField( fphys2d(:,:,3) - fm(:,:,2) ) );
end