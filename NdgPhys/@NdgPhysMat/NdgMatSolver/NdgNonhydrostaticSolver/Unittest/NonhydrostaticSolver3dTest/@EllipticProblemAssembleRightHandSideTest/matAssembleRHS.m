function RHS = matAssembleRHS( obj, fphys, deltatime, PUPX, PUPS, PVPY, PVPS, PWPS, PSPX, PSPY )
%ASSEMBLERHS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
mesh = obj.meshUnion;
cell = mesh.cell;
TempRHS = obj.NonhydrostaticSolver.rho/deltatime*(PUPX + PUPS.*PSPX + PVPY + PVPS.*PSPY + 1./fphys{1}(:,:,4).*PWPS);
RHS = zeros(cell.Np, mesh.K);
for i = 1:mesh.K
    RHS(:,i) = diag(mesh.J(:,i))*cell.M*TempRHS(:,i);
end
end

