function matAssembleGlobalStiffMatrix(obj)
x = obj.NonhydrostaticSolver.BoundaryEdge.xb;
DirichletData = eval(obj.Cexact);

TempRHS = obj.RHS;
Tau = 30;

[ obj.RHS, obj.StiffMatrix ] = mxAssembleGlobalStiffMatrixWithBCsImposed(obj.NonhydrostaticSolver.SPNPX, ...
   obj.RHS, DirichletData, obj.NonhydrostaticSolver.BoundaryEdge2d, obj.NonhydrostaticSolver.BoundaryEdge,...
   obj.NonhydrostaticSolver.cell, obj.NonhydrostaticSolver.mesh, int8(obj.meshUnion.mesh2d.BoundaryEdge.ftype) );

mesh2d = obj.meshUnion.mesh2d;
mesh3d = obj.meshUnion;
cell = obj.meshUnion.cell;
nx = obj.NonhydrostaticSolver.BoundaryEdge.nx(1);
ny = obj.NonhydrostaticSolver.BoundaryEdge.ny(1);

Dx = diag(mesh3d.rx(:,1))*cell.Dr + diag(mesh3d.sx(:,1))*cell.Ds;
Dy = diag(mesh3d.ry(:,1))*cell.Dr + diag(mesh3d.sy(:,1))*cell.Ds;

M2d = diag(obj.NonhydrostaticSolver.BoundaryEdge.Js(:,1))*mesh2d.cell.M;
WeightedM2d = M2d*diag(DirichletData(:,1));

Index = cell.Fmask(:,1);
Contribution1 = nx .* ( Dx(Index,:)'*WeightedM2d ) + ny .* ( Dy(Index,:)'*WeightedM2d );

Contribution = sum(Contribution1,2);

Contribution(Index) = Contribution(Index) - sum(Tau*WeightedM2d, 2);
InvM = inv(diag(mesh3d.J(:,1))*cell.M);
TempRHS(1:cell.Np) = TempRHS(1:cell.Np) + (InvM * Contribution)';

end