clear,clc;
disp('************************For partial derivative calculation************************')
Solver = EllipticProblemCalculatePartialDerivative3d(1,1,50,5);
Solver.EllipticProblemSolve;
clear;
disp('************************End partial derivative calculation************************')

disp('************************For Global Stiff Matrix Assemble************************')
Solver = EllipticProblemMatrixAssembleTest3dNew(1,1,10,5);
Solver.EllipticProblemSolve;
clear;
disp('************************End Global Stiff Matrix Assemble************************')

disp('************************For Right Hand Side Calculation************************')
Solver = EllipticProblemAssembleRightHandSideTest(1,1,50,5);
Solver.EllipticProblemSolve;
clear;
disp('************************End Right Hand Side Calculation************************')

disp('************************For impose boundary condition************************')
Solver = NonhydroImposeBCsTest3d(1, 2, 10, 5);
Solver.EllipticProblemSolve;
clear;
disp('************************End impose boundary condition************************')

disp('************************For assemble the final global stiff matrix************************')
Solver = EllipticProblemAssembleFinalStiffMatrix3d(1,1,10,3);
Solver.EllipticProblemSolve;
clear;
disp('************************End assemble the final global stiff matrix************************')

disp('************************For update conservative final velocity************************')
Solver = EllipticProblemUpdateConservativeFinalVelocityTest(1,1,50,5);
Solver.EllipticProblemSolve;
clear;
disp('************************End update conservative final velocity************************')

