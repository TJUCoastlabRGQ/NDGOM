function matCalculateFphysDerivative(obj, mesh, fphys, physClass)



[ obj.fhx ] = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(fphys{1}(:,:,1) ,[1 2]));

obj.HBx = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(fphys{1}(:,:,1) + 2*fphys{1}(:,:,4)  ,[1 2]));

% obj.fhx + 2 * physClass.zGrad{1}(:,:,1);

[ obj.hux ] = obj.matCalculateUpwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(fphys{1}(:,:,2) ,[1 2]));

[ obj.H2Bx ] = obj.matCalculateDownwindedFphysDerivative( mesh, fphys, physClass, ...
    num2cell(obj.HBx ,[1 2]));

end
