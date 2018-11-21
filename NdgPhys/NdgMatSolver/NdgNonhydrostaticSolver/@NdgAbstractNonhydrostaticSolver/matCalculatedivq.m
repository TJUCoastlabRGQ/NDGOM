function divq = matCalculatedivq(obj, mesh, Qx, Qy)
Qxr = mesh.cell.Dr * Qx; Qxs = mesh.cell.Ds * Qx;
Qyr = mesh.cell.Dr * Qy; Qys = mesh.cell.Ds * Qy;
divq = mesh.rx.*Qxr + mesh.sx .* Qxs + mesh.ry.*Qyr + mesh.sy .* Qys;
end