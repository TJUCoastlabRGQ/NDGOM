function matEvaluateSourceTerm( obj, fphys, fphys2d, time )
% frhs = frhs + BottomTerm
obj.matEvaluateTopographySourceTerm( fphys, fphys2d, time );
% frhs = frhs + CoriolisTerm
obj.coriolisSolver.evaluateCoriolisTermRHS(obj, fphys);
end