function penaltyterm = matEvaluatePenaltyTerm(obj, mesh, fm, fp, directionalVector)
hmin = 2*mesh.J(mesh.eidM)./mesh.Js;
tau = mesh.cell.Np./hmin;
penaltyterm = fm .* directionalVector .* tau -  fp .* directionalVector .* tau;
end