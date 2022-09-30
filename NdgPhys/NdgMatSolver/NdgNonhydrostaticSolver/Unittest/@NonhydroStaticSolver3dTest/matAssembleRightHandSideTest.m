function matAssembleRightHandSideTest( obj )

CversionRHS = obj.NonhydrostaticSolver.TestAssembleRHS( obj, obj.fphys, obj.fphys2d, 0.5 );

CversionRHS = full( CversionRHS );

matVersionRHS = 1000/0.5*(obj.NonhydrostaticSolver.PUPX + obj.NonhydrostaticSolver.PUPS.*obj.NonhydrostaticSolver.PSPX + ...
    obj.NonhydrostaticSolver.PVPY + obj.NonhydrostaticSolver.PVPS.*obj.NonhydrostaticSolver.PSPY + ...
    1./obj.fphys{1}(:,:,4).*obj.NonhydrostaticSolver.PWPS);

disp("=====================For Assembling right hand side part==============================");
disp("The maximum difference is:");
disp(max(max(CversionRHS(:) -  matVersionRHS(:))));
disp("The minimum difference is:");
disp(min(min(CversionRHS(:) -  matVersionRHS(:))));
disp("The maximum value of the C version is:");
disp(max(max(CversionRHS(:))));
disp("The minimum value of the C version is:");
disp(min(min(CversionRHS(:))));
disp("======================================================================================");
end