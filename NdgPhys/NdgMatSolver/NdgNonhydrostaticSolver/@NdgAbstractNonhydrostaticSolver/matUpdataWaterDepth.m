function WaterDepth  = matUpdataWaterDepth(obj, physClass, mesh, fphys, deltafluxhu, deltafluxhv)
deltafluxhuP = deltafluxhu(mesh.eidP);deltafluxhvP = deltafluxhv(mesh.eidP);
index = ( mesh.eidtype ~= 0 );
deltafluxhuP(index) = - deltafluxhuP( index );deltafluxhvP(index) = deltafluxhvP( index );
fluxq = mesh.nx .* (deltafluxhu(mesh.eidM) - deltafluxhuP)/2 + mesh.ny .* (deltafluxhv(mesh.eidM) - deltafluxhvP)/2;
WaterDepth =  - physClass.dt * (mesh.rx .* (mesh.cell.Dr * deltafluxhu) + mesh.sx .* (mesh.cell.Ds * deltafluxhu) + ...
   mesh.ry .* (mesh.cell.Dr * deltafluxhv) + mesh.sy .* (mesh.cell.Ds * deltafluxhv) - mesh.cell.LIFT*(mesh.Js .* fluxq)./mesh.J) + fphys{1}(:,:,1);
end