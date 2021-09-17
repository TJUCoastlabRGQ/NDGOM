function [ FinalStiffMatrix ] = matAssembleFinalGlobalStiffMatrix( obj, phpx, phpy )
%MATASSEMBLEFINALGLOBALSTIFFMATRIX 此处显示有关此函数的摘要
%   此处显示详细说明
FinalStiffMatrix = full(obj.StiffMatrix);
PNPX = full(obj.NonhydrostaticSolver.PNPX);
PNPY = full(obj.NonhydrostaticSolver.PNPY);
PNPS = full(obj.NonhydrostaticSolver.PNPS);
zx = obj.fphys{1}(:,:,8);
zy = obj.fphys{1}(:,:,9);

Depth = obj.fphys{1}(:,:,4);

StiffTerm = bsxfun(@times, 1./Depth(:).*phpx(:), PNPX) + bsxfun(@times, 1./Depth(:).*phpy(:), PNPY) - ...
    bsxfun(@times, (1./Depth(:).*zx(:) + (1 + obj.meshUnion.z(:))./Depth(:).*phpx(:)).*(1./Depth(:).*phpx(:)), PNPS) - ...
    bsxfun(@times, (1./Depth(:).*zy(:) + (1 + obj.meshUnion.z(:))./Depth(:).*phpy(:)).*(1./Depth(:).*phpy(:)), PNPS);

Np = obj.meshUnion.cell.Np;
MassMatrix = diag(obj.meshUnion.J(:,1))*obj.meshUnion.cell.M;
for i = 1:obj.meshUnion.K * Np
    for j = 1:obj.meshUnion.K
        TempData = StiffTerm((i-1)*obj.meshUnion.K*Np + (j-1)*Np + 1:(i-1)*obj.meshUnion.K*Np + j*Np);
        FinalStiffMatrix((i-1)*obj.meshUnion.K*Np + (j-1)*Np + 1:(i-1)*obj.meshUnion.K*Np + j*Np) = ...
            FinalStiffMatrix((i-1)*obj.meshUnion.K*Np + (j-1)*Np + 1:(i-1)*obj.meshUnion.K*Np + j*Np) + (MassMatrix * TempData')';
    end
end

end

