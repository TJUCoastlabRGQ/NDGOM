function matAverageBoundaryTopoGrad(obj, PhysClass, mesh)
%> @brief Function to set the averaged bottom topography gradient
%> @details Function to set the averaged bottom topography gradient
%> @param[in] physClass The hydrostatic solver
%> @param[in] mesh The mesh object
tempbx = PhysClass.zGrad{1}(:,:,1);tempby = PhysClass.zGrad{1}(:,:,2);
obj.bx = tempbx; obj.by = tempby;
% obj.bx(mesh.eidM) = (tempbx(mesh.eidM) + tempbx(mesh.eidP))/2;
% obj.by(mesh.eidM) = (tempby(mesh.eidM) + tempby(mesh.eidP))/2;
end