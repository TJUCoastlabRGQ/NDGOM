function coriolisSolver3d = initcoriolisSolver3d(obj)
if obj.option.isKey('CoriolisType') % the option exist
     %doing nothing
else % the option does not exist
    coriolisSolver3d = NonCoriolisTermSolver();
end
end