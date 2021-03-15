clear, clc;

%% This is test part for term $\frac{\partial^2 p}{\partial \sigma^2}$
% Solver = EllipticProblem3d( 1, 1, 1, 100 );
% matSolver = EllipticProblemTest3d(1, 1, 2, 100);
% matSolver.EllipticProblemSolve;
syms x y z;
% Exact = sin(-pi/2*x);
% FDiffExact = diff(Exact, x);
% DiffExact = diff(diff(Exact, x),x);
% x = Solver.meshUnion.z;
% ExactValue = eval(Exact);
% data = eval(DiffExact);
% 
% Re = Solver.NonhydrostaticSolver.SPNPS\data(:);
% Residual = Re - ExactValue(:);
% clear matSolver;
% matSolver = FirstOrderProblemTest3d(1, 1, 2, 100);
% data = eval(FDiffExact);
% Re = Solver.NonhydrostaticSolver.PNPS\data(:);
% Residual = Re - ExactValue(:);
% 
Solver = PartailDeriveTest;
uE = x*(20-x)*y*(20-y)*z;
vE = x*(20-x)*y*(20-y)*z;
hE = x*(20-x)*y*(20-y) + 10;
etaE = x*(20-x)*y*(20-y);

PSPX = -1*z/hE*diff(hE,x)-1/hE*diff(etaE,x);
PSPY = -1*z/hE*diff(hE,y)-1/hE*diff(etaE,y);
SPSPX = diff(PSPX, x);
SPSPY = diff(PSPY, y);
% SPSPX = z/hE/hE*diff(hE,x)*diff(hE,x) - ...
%     z/hE*diff(diff(hE,x),x)+1/hE/hE*diff(hE,x) * diff(etaE,x) - ...
%     1/hE*diff(diff(etaE,x),x);
% SPSPY = z/hE/hE*diff(hE,y)*diff(hE,y) - ...
%     z/hE*diff(diff(hE,y),y)+1/hE/hE*diff(hE,y) * diff(etaE,y) - ...
%     1/hE*diff(diff(etaE,y),y);
PUPX = diff(uE,x);
PVPY = diff(vE,y);
PUPS = diff(uE,z);
PVPS = diff(vE,z);

x = Solver.meshUnion.x;
y = Solver.meshUnion.y;
z = Solver.meshUnion.z;
fphys = Solver.fphys;
fphys2d = Solver.fphys2d;
fphys{1}(:,:,1) = eval(hE).*eval(uE);
fphys{1}(:,:,2) = eval(hE).*eval(vE);
Solver.NonhydrostaticSolver.TestPartialDerivativeCalculation( Solver, fphys, fphys2d, 1);
PSPXE = eval(PSPX);
PSPYE = eval(PSPY);
SPSPXE = eval(SPSPX);
SPSPYE = eval(SPSPY);
PUPXE = eval(PUPX);
PVPYE = eval(PVPY);
PUPSE = eval(PUPS);
PVPSE = eval(PVPS);

r1 = Solver.NonhydrostaticSolver.PSPX - PSPXE;
r2 = Solver.NonhydrostaticSolver.PSPY - PSPYE;
r3 = Solver.NonhydrostaticSolver.SPSPX - SPSPXE; 
r4 = Solver.NonhydrostaticSolver.SPSPY - SPSPYE; 
r5 = Solver.NonhydrostaticSolver.PUPX - PUPXE;
r6 = Solver.NonhydrostaticSolver.PVPY - PVPYE;
r7 = Solver.NonhydrostaticSolver.PUPS - PUPSE;
r8 = Solver.NonhydrostaticSolver.PVPS - PVPSE;


%% This is to assemble data into the sparse matrix, checked
disp('======================================================================');
InsertData = [1 2 3 4; 5 6 7 8; 1 4 5 6; 7 8 9 6];
InsertRow = [1 3 4 6];
Sdata = zeros(7,4);
Sdata(InsertRow,:) = InsertData;
InputSparse = sparse(Sdata);
OutputSparse = AssembleDataIntoSparseMatrix(7,4,InsertData, InputSparse);
disp(InputSparse - OutputSparse);
disp('Assemble data into sparse matrix has been verified.');
%% This is to fetch the data from sparse matrix, checked
% disp('======================================================================');
% InsertData = [1 2 3 4; 5 6 7 8; 1 4 5 6; 7 8 9 6];
% InsertRow = [1 3 4 6];
% Sdata = zeros(7,4);
% Sdata(InsertRow,:) = InsertData;
% InputSparse = sparse(Sdata);
% OutputData = FetchDataFromSparseMatrix(InputSparse,4,4);
% disp(InsertData - OutputData);
% disp('Fetch data from sparse matrix has been verified.');

%% This is to fetch the diff matrix, checked
disp('======================================================================');
Solver = StandingWaveInAClosedChannel( 1, 1, 50, 2 );
Ele = 47;
fetchRow = Solver.meshUnion.cell.Fmask(:,end);
Dr = Solver.meshUnion.cell.Dr;
Ds = Solver.meshUnion.cell.Ds;
Dt = Solver.meshUnion.cell.Dt;
Np3d = Solver.meshUnion.cell.Np;
M = Solver.meshUnion.cell.M;
M2d = Solver.meshUnion.mesh2d.cell.M;
rx = Solver.meshUnion.rx;
sx = Solver.meshUnion.sx;
ry = Solver.meshUnion.ry;
sy = Solver.meshUnion.sy;
tz = Solver.meshUnion.tz;
J = Solver.meshUnion.J;
J2d = Solver.meshUnion.mesh2d.J;
EDx = diag(rx(:,Ele))*Dr + diag(sx(:,Ele))*Ds;
EDy = diag(ry(:,Ele))*Dr + diag(sy(:,Ele))*Ds;
EDz = diag(tz(:,Ele))*Dt;
EFDx2d = EDx(fetchRow,:);
EFDy2d = EDy(fetchRow,:);
EM3d = diag(J(:,Ele))*M;
EInvM3d = inv(EM3d);
EDzTM3d = EDz'*EM3d;
EMIMDzTM3dDz = -1*EInvM3d*EDz'*EM3d*EDz;
EM2d = diag(J2d(:,Ele))*M2d;
EFacialDiffxMultM2d = EFDx2d'*EM2d;
EFacialDiffyMultM2d = EFDy2d'*EM2d;
EM2dMultFacialDiffx = EM2d*EFDx2d;
EM2dMultFacialDiffy = EM2d*EFDy2d;
ETestTempContributionx = zeros(Np3d);
ETestTempContributionx(fetchRow,:) = EM2dMultFacialDiffx;
ETestTempContributiony = zeros(Np3d);
ETestTempContributiony(fetchRow,:) = EM2dMultFacialDiffy;
EDirichletBuff2 = zeros(Np3d);
EDirichletBuff2(fetchRow,fetchRow) = -1*EM2d;
TempContributionBuff2 = zeros(Np3d);
TempContributionBuff2(:,fetchRow) = 0.5*EFDx2d'*EM2d;
% OutputDiffMatrix = FetchDataFromDiffMatrix(fetchRow, Dr);
warning('off');
[Dx, Dy, Dz, FDx2d, FDy2d, M3d, InvM3d, DzTM3d, MIMDzTM3dDz, M2d, FacialDiffxMultM2d, ...
    FacialDiffyMultM2d, M2dMultFacialDiffx, M2dMultFacialDiffy, TestTempContributionx, ...
    TestTempContributiony, DirichletBuff2, ETempContributionBuff2] = FetchDataFromDiffMatrix(fetchRow, struct(Solver.meshUnion), ...
    struct(Solver.meshUnion.cell), struct(Solver.meshUnion.mesh2d), struct(Solver.meshUnion.mesh2d.cell), Ele);
warning('on');
disp(Dx - EDx);
disp('Diff matrix in x direction has been verified.');
disp(Dy - EDy);
disp('Diff matrix in y direction has been verified.');
disp(Dz - EDz);
disp('Diff matrix in z direction has been verified.');
disp(FDx2d - EFDx2d);
disp('Diff matrix for facial diffmatrix in x direction has been verified.');
disp(FDy2d - EFDy2d);
disp('Diff matrix for facial diffmatrix in y direction has been verified.');
disp(M3d - EM3d);
disp('Three dimensional mass matrix has been verified.');
disp(InvM3d - EInvM3d);
disp('Three dimensional inverse mass matrix has been verified.');
disp(DzTM3d - EDzTM3d);
disp('Multiplication of transpose diff matrix with mass matrix has been verified.');
disp(MIMDzTM3dDz - EMIMDzTM3dDz);
disp('Multiplication of several matrix has been verified.');
disp(EM2d - M2d);
disp('Two dimensional mass matrix has been verified.');
disp(EFacialDiffxMultM2d - FacialDiffxMultM2d);
disp('Multiplication of two-dimensional diff matrix in x direction with two-dimensional mass matrix has been verified.');
disp(EFacialDiffyMultM2d - FacialDiffyMultM2d);
disp('Multiplication of two-dimensional diff matrix in y direction with two-dimensional mass matrix has been verified.');
disp(EM2dMultFacialDiffx - M2dMultFacialDiffx);
disp('Multiplication of two-dimensional mass matrix in x direction with two-dimensional diff matrix has been verified.');
disp(EM2dMultFacialDiffy - M2dMultFacialDiffy);
disp('Multiplication of two-dimensional mass matrix in y direction with two-dimensional diff matrix has been verified.');
disp(ETestTempContributionx - TestTempContributionx);
disp('Assemble contribution into row in x direction has been verified.');
disp(ETestTempContributiony - TestTempContributiony);
disp('Assemble contribution into row in y direction has been verified.');
disp(EDirichletBuff2 - DirichletBuff2);
disp('Assemble contribution into row and column has been verified.');
disp(ETempContributionBuff2 - TempContributionBuff2);
disp('Assemble contribution into column has been verified.');
