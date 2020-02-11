function RHS = matEvaluateConservativeNonhydrostaticRHS(obj, fphys, PhysClass)

RHS = 2 * fphys{1}(:,:,6) + fphys{1}(:,:,1) .* obj.hux + fphys{1}(:,:,1) .* obj.hvy - fphys{1}(:,:,2) .* obj.HBx ...
    - fphys{1}(:,:,3) .* obj.HBy;

end