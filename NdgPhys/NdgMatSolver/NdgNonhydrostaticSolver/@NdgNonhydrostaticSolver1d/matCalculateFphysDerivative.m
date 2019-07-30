function matCalculateFphysDerivative(obj, mesh, fphys, physClass)

[ obj.HBx ] = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...    % abstract in nonhydrostatic solver 2d and 1d, different for Gauss quad and quad free version
    num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,3),[1 2]));

[ obj.fhx ] = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(fphys{1}(:,:,1) ,[1 2]));

[ obj.hux ] = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(fphys{1}(:,:,2) ,[1 2]));

[ obj.H2Bx ] = obj.matCalculateDownwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(obj.HBx ,[1 2]) );

end
