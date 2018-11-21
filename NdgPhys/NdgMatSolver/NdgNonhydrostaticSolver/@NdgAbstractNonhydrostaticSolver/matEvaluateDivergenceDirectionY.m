function divY = matEvaluateDivergenceDirectionY(obj, mesh, Qy)
divY = mesh.ry .* (mesh.cell.Dr * Qy) + mesh.sy .* (mesh.cell.Ds * Qy);
end