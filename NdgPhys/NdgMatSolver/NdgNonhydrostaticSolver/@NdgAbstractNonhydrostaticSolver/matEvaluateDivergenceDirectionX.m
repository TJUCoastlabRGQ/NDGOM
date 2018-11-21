function divX = matEvaluateDivergenceDirectionX(obj, mesh, Qx)
divX = mesh.rx .* (mesh.cell.Dr * Qx) + mesh.sx .* (mesh.cell.Ds * Qx);
end