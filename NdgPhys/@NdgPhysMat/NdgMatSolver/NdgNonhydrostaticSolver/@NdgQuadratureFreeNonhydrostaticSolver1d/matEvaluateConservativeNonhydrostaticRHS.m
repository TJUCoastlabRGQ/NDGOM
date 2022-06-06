function RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, PhysClass)

RHS = 2 * fphys{1}(:,:,5) + fphys{1}(:,:,1) .* obj.hux  - fphys{1}(:,:,2) .* obj.HBx;

end