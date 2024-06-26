function matCalculateFphysDerivative(obj, mesh, fphys, physClass)

[obj.HBx, obj.HBy] = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...    % abstract in nonhydrostatic solver 2d and 1d, different for Gauss quad and quad free version
    num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), num2cell(fphys{1}(:,:,1) +...
    2 * fphys{1}(:,:,4),[1 2]));

[obj.fhx, obj.fhy] = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(fphys{1}(:,:,1) ,[1 2]), num2cell(fphys{1}(:,:,1),[1 2]));

[obj.hux, obj.hvy] = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(fphys{1}(:,:,2) ,[1 2]), num2cell(fphys{1}(:,:,3),[1 2]));

[obj.H2Bx, obj.H2By] = obj.matCalculateDownwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(obj.HBx ,[1 2]), num2cell(obj.HBy ,[1 2]));

end
