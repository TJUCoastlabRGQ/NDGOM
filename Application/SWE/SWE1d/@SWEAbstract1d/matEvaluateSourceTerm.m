function matEvaluateSourceTerm( obj, fphys )
% frhs = frhs + BottomTerm
obj.matEvaluateTopographySourceTerm( fphys );

% frhs = frhs + FrictionTerm
obj.frictionSolver.evaluateFrictionTermRHS(obj, fphys);

obj.NonhydrostaticSolver.evaluateNonhydroRHS(obj, fphys);
end