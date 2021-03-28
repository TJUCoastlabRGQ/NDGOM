function matAssembleGlobalStiffMatrix(obj)
% x = obj.NonhydrostaticSolver.BoundaryEdge.xb;
% y = obj.NonhydrostaticSolver.BoundaryEdge.yb;
% DirichletData = eval(obj.Cexact);
% 
% [ obj.RHS, obj.StiffMatrix ] = mxAssembleGlobalStiffMatrixWithBCsImposed(obj.NonhydrostaticSolver.SPNPX + obj.NonhydrostaticSolver.SPNPY, ...
%    obj.RHS, DirichletData, obj.NonhydrostaticSolver.BoundaryEdge2d, obj.NonhydrostaticSolver.BoundaryEdge,...
%    obj.NonhydrostaticSolver.cell, obj.NonhydrostaticSolver.mesh, int8(obj.meshUnion.mesh2d.BoundaryEdge.ftype) );
obj.StiffMatrix = obj.NonhydrostaticSolver.SPNPS;
end