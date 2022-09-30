function matAssembleGlobalStiffMatrixTest( obj )

MatVersionStiffMatrix = zeros( obj.meshUnion.K * obj.meshUnion.cell.Np );
SPNPX = full(obj.NonhydrostaticSolver.SPNPX);
SPNPY = full(obj.NonhydrostaticSolver.SPNPY);
PNPS  = full(obj.NonhydrostaticSolver.PNPS);
SPNPS  = full(obj.NonhydrostaticSolver.SPNPS);
MSPNPX = full(obj.NonhydrostaticSolver.MSPNPX);
MSPNPY = full(obj.NonhydrostaticSolver.MSPNPY);
InvHeightSquare = 1./obj.fphys{1}(:,:,4)./obj.fphys{1}(:,:,4);
for i = 1:obj.meshUnion.K * obj.meshUnion.cell.Np
    MatVersionStiffMatrix(:,i) = SPNPX(:, i) + SPNPY(:, i) + ...
        ( obj.NonhydrostaticSolver.SQPSPX(:) + obj.NonhydrostaticSolver.SQPSPY(:) + InvHeightSquare(:) ).*SPNPS(:,i) + ...
        2 * ( obj.NonhydrostaticSolver.PSPX(:).*MSPNPX(:,i) + obj.NonhydrostaticSolver.PSPY(:).*MSPNPY(:,i) ) + ...
        ( obj.NonhydrostaticSolver.SPSPX(:) + obj.NonhydrostaticSolver.SPSPY(:) ) .* PNPS(:,i);
end

CVersionStiffMatrix = obj.NonhydrostaticSolver.TestAssembleGlobalStiffMatrix( obj, obj.fphys, obj.fphys2d);

CVersionStiffMatrix = full(CVersionStiffMatrix);

disp("==========================For Stiff Matrix==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(MatVersionStiffMatrix-CVersionStiffMatrix)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(MatVersionStiffMatrix-CVersionStiffMatrix)));
disp("The maximum value of stiff matrix of the C version is:");
disp(max(max(CVersionStiffMatrix)));
disp("The minimum value of stiff matrix of the C version is:");
disp(min(min(CVersionStiffMatrix)));
disp("=============================================================================");

end