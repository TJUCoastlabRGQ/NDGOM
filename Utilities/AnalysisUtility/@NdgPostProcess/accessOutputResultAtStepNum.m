function [ fphys ] = accessOutputResultAtStepNum(obj, stepId)

fphys = cell( obj.Nmesh, 1 );
for m = 1:obj.Nmesh
    Np = obj.meshUnion(m).mesh2d.cell.Np;
    K = obj.meshUnion(m).mesh2d.K;
    fphys{m} = ncread( obj.outputFile{m}, 'fphys', [1, 1, 1, stepId], [Np, K, obj.Nvar, 1]);
%     fphys{m} = ncread( obj.outputFile{m}, 'fphys2d', [1, 1, 1, stepId], [Np, K, 3, 1]);
end

end% func