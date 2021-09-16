function matAssembleStiffMatrixAndRHS( obj, TempStiffMatrix, TempRHS, PSPX, PSPY, Wold, Wnew, deltatime, PWPS, Unew, Uold, Vnew, Vold, PUPX, PUPY, PUPS, PVPX, PVPY, PVPS, PHPX, PHPY)

rho = obj.NonhydrostaticSolver.rho;
K31 = PSPX;
K32 = PSPY;

mesh = obj.meshUnion(1).mesh2d;
cell2d = mesh.cell;
edge = mesh.InnerEdge;
W = cell(1);
W{1} = Wnew;
% For term $\frac{\partial w}{\partial x}$
[ fm, fp ] = edge.matEvaluateSurfValue( W );
fluxMx = edge.nx .* fm;
fluxPx = edge.nx .* fp;
fluxSx = edge.nx .* (fm + fp)./2;
rhsX = edge.matEvaluateStrongFormEdgeRHS( fluxMx, fluxPx, fluxSx );
Wx = mesh.rx .* (cell2d.Dr * Wnew) + mesh.sx .* (cell2d.Ds * Wnew) - rhsX;

obj.MatWx = Wx;

% For term $\frac{\partial w}{\partial x}$
fluxMy = edge.ny .* fm;
fluxPy = edge.ny .* fp;
fluxSy = edge.ny .* (fm + fp)./2;
rhsY = edge.matEvaluateStrongFormEdgeRHS( fluxMy, fluxPy, fluxSy );
Wy = mesh.ry .* (cell2d.Dr * Wnew) + mesh.sy .* (cell2d.Ds * Wnew) - rhsY;

obj.MatWy = Wy;

edge3d = obj.meshUnion.BottomBoundaryEdge;
CPWPS = cell(1);
CPWPS{1} = PWPS;
[ BotBEPWPS, ~ ] = edge3d.matEvaluateSurfValue( CPWPS );

CPUPX = cell(1);
CPUPX{1} = PUPX;
[ BotBEPUPX, ~ ] = edge3d.matEvaluateSurfValue( CPUPX );

CPUPY = cell(1);
CPUPY{1} = PUPY;
[ BotBEPUPY, ~ ] = edge3d.matEvaluateSurfValue( CPUPY );

CPUPS = cell(1);
CPUPS{1} = PUPS;
[ BotBEPUPS, ~ ] = edge3d.matEvaluateSurfValue( CPUPS );

CPVPX = cell(1);
CPVPX{1} = PVPX;
[ BotBEPVPX, ~ ] = edge3d.matEvaluateSurfValue( CPVPX );

CPVPY = cell(1);
CPVPY{1} = PVPY;
[ BotBEPVPY, ~ ] = edge3d.matEvaluateSurfValue( CPVPY );

CPVPS = cell(1);
CPVPS{1} = PVPS;
[ BotBEPVPS, ~ ] = edge3d.matEvaluateSurfValue( CPVPS );


CPHPX = cell(1);
CPHPX{1} = PHPX;
[ BotBEPHPX, ~ ] = edge3d.matEvaluateSurfValue( CPHPX );

CPHPY = cell(1);
CPHPY{1} = PHPY;
[ BotBEPHPY, ~ ] = edge3d.matEvaluateSurfValue( CPHPY );

[fm, ~] = edge3d.matEvaluateSurfValue(obj.fphys);
WaterDepth = fm(:,:,4);
InverseHeight = 1./WaterDepth;
pzpx = fm(:,:,8);
pzpy = fm(:,:,9);

K33 = K31.^2 + K32.^2 + 1./obj.fphys{1}(:,:,4)./obj.fphys{1}(:,:,4);

Coe = cell(1);
Coe{1}(:,:,1) = K31;
Coe{1}(:,:,2) = K32;
Coe{1}(:,:,3) = K33;

[BotBECoe, ~] = edge3d.matEvaluateSurfValue(Coe);

PNPS = -1*rho.*WaterDepth.*( (Wnew-Wold)./deltatime+Unew.*(Wx - InverseHeight.*BotBEPWPS.*pzpx) + Vnew.*(Wy - InverseHeight.*BotBEPWPS.*pzpy) + Wnew .* InverseHeight .* BotBEPWPS );

PNPX = -1*rho.*( (Unew-Uold)./deltatime+Unew.*(BotBEPUPX - InverseHeight.*BotBEPUPS.*pzpx) + Vnew.*(BotBEPUPY - InverseHeight.*BotBEPUPS.*pzpy) + Wnew .* InverseHeight .* BotBEPUPS + obj.gra*BotBEPHPX + InverseHeight.*pzpx.*(obj.gra*WaterDepth-1/rho.*PNPS));

PNPY = -1*rho.*( (Vnew-Vold)./deltatime+Unew.*(BotBEPVPX - InverseHeight.*BotBEPVPS.*pzpx) + Vnew.*(BotBEPVPY - InverseHeight.*BotBEPVPS.*pzpy) + Wnew .* InverseHeight .* BotBEPVPS + obj.gra*BotBEPHPY + InverseHeight.*pzpy.*(obj.gra*WaterDepth-1/rho.*PNPS));

obj.MatPNPX = (PNPX .* BotBECoe(:,:,1));

obj.MatPNPY = (PNPY .* BotBECoe(:,:,2));

obj.MatPNPS = (PNPS .* BotBECoe(:,:,3));

PNPX = diag(edge3d.Js(:,1))*edge3d.M * (PNPX .* BotBECoe(:,:,1));

PNPY = diag(edge3d.Js(:,1))*edge3d.M * (PNPY .* BotBECoe(:,:,2));

PNPS = diag(edge3d.Js(:,1))*edge3d.M * (PNPS .* BotBECoe(:,:,3));

Index = bsxfun(@plus,(edge3d.FToE(1,:) - 1)*obj.meshUnion.cell.Np, edge3d.FToN1);

obj.RHS = TempRHS;

obj.RHS(Index) = TempRHS(Index) - (-1*(PNPS + PNPY + PNPX));

end