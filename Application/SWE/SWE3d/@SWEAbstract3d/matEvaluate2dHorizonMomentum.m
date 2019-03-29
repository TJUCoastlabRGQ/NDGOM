function fphys2d = matEvaluate2dHorizonMomentum( obj, mesh3d, fphys2d, fphys3d )
%> @brief Function used to calculate the vertically averaged horizontal momentum term
%> 
%> More detailed description.
%>
%> @param mesh3d The three dimensional mesh object
%> @param fphys2d The two dimensional physical field
%> @param fphys3d The three dimensional physical field
%>
%> @fphys2d fphys2d The two dimensional physical field with the vertically averaged horizontal momentum term updated


fphys2d(:, :, 2) = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 1) );
fphys2d(:, :, 3) = mesh3d.VerticalColumnIntegralField( fphys3d(:, :, 2) );

end