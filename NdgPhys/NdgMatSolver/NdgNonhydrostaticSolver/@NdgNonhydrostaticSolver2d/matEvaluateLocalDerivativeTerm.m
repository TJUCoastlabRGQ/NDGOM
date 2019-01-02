function [TermX, TermY] = matEvaluateLocalDerivativeTerm(obj, mesh, fieldValue)
TermX = mesh.rx .* (mesh.cell.Dr * fieldValue) + mesh.sx .* (mesh.cell.Ds * fieldValue);
TermY = mesh.ry .* (mesh.cell.Dr * fieldValue) + mesh.sy .* (mesh.cell.Ds * fieldValue);
end