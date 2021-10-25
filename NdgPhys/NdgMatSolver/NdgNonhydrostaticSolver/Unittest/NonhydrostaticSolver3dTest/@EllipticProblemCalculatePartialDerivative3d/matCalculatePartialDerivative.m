function matCalculatePartialDerivative( obj )
%MATCALCULATEPARTIALDERIVATIVE 此处显示有关此函数的摘要
%   此处显示详细说明
InnerEdge = obj.meshUnion.InnerEdge;
BoundaryEdge = obj.meshUnion.BoundaryEdge;
BottomEdge = obj.meshUnion.BottomEdge;
BottomBoundaryEdge = obj.meshUnion.BottomBoundaryEdge;
SurfaceBoundaryEdge = obj.meshUnion.SurfaceBoundaryEdge;

[obj.matPUPX, obj.matPUPY, obj.matPVPX, obj.matPVPY, obj.matPHPX, obj.matPHPY] ...
    = CalculateHorizontalDeriveTerm( obj, InnerEdge, BoundaryEdge);

[obj.matPUPS, obj.matPVPS, obj.matPWPS, obj.matWnew, obj.matUnew, obj.matVnew ] = CalculateVerticalDerivativeTerm( obj, BottomEdge, BottomBoundaryEdge, SurfaceBoundaryEdge );

[obj.matPSPX, obj.matPSPY, obj.matSQPSPX, obj.matSQPSPY] = CalculateSigmaRelatedTerm( obj );

end

function [ PSPX, PSPY, SQPSPX, SQPSPY ] = CalculateSigmaRelatedTerm( obj )
PSPX = -1 * 1./obj.fphys{1}(:,:,4) .* obj.fphys{1}(:,:,8) - (1 + obj.meshUnion.z)./obj.fphys{1}(:,:,4) .* obj.matPHPX;
PSPY = -1 * 1./obj.fphys{1}(:,:,4) .* obj.fphys{1}(:,:,9) - (1 + obj.meshUnion.z)./obj.fphys{1}(:,:,4) .* obj.matPHPY;
SQPSPX = PSPX.^2;
SQPSPY = PSPY.^2;
end

function [PUPS, PVPS, PWPS, Wnew, Unew, Vnew ] = CalculateVerticalDerivativeTerm( obj, BottomEdge, BottomBoundaryEdge, SurfaceBoundaryEdge )
mesh = obj.meshUnion;
Cell = mesh.cell;
[Surffm, ~] = SurfaceBoundaryEdge.matEvaluateSurfValue( obj.fphys );
ueta = Surffm(:,:,1)./Surffm(:,:,4);
veta = Surffm(:,:,2)./Surffm(:,:,4);
TempData = cell(1);
TempData{1}(:,:,1) = obj.matPHPX;
TempData{1}(:,:,2) = obj.matPHPY;
[petapdfm, ~] = SurfaceBoundaryEdge.matEvaluateSurfValue( TempData );
% $w_{\eta} = \frac{\partial \eta}{\partial t} + u_{\eta}\frac{\partial \eta}{\partial x} + v_{\eta}\frac{\partial \eta}{\partial y}$
%$\frac{\partial \eta}{\partial x} = \frac{\partial H}{\partial x} + \frac{\partial z_b}{\partial x}$
Weta = obj.frhs2d{1} + ueta .* (petapdfm(:,:,1) + Surffm(:,:,8)) + veta .* (petapdfm(:,:,2) + Surffm(:,:,9));

Hu = obj.fphys{1}(:,:,1);
Hv = obj.fphys{1}(:,:,2);
Hw = obj.fphys{1}(:,:,11);
H = obj.fphys{1}(:,:,4);
u = Hu./H;
v = Hv./H;
w = Hw./H;

PUPS = mesh.tz .* (Cell.Dt * u);
PVPS = mesh.tz .* (Cell.Dt * v);
PWPS = mesh.tz .* (Cell.Dt * w);

[fm, fp] = BottomEdge.matEvaluateSurfValue( obj.fphys );
PUPSBotEFluxM = BottomEdge.nz .* fm(:,:,1) ./ fm(:,:,4);
PUPSBotEFluxP = BottomEdge.nz .* fp(:,:,1) ./ fp(:,:,4);
PUPSBotEFluxS = ( PUPSBotEFluxM + PUPSBotEFluxP ) ./ 2;
PUPS = PUPS - BottomEdge.matEvaluateStrongFormEdgeRHS( PUPSBotEFluxM, PUPSBotEFluxP, PUPSBotEFluxS );

PVPSBotEFluxM = BottomEdge.nz .* fm(:,:,2) ./ fm(:,:,4);
PVPSBotEFluxP = BottomEdge.nz .* fp(:,:,2) ./ fp(:,:,4);
PVPSBotEFluxS = ( PVPSBotEFluxM + PVPSBotEFluxP ) ./ 2;
PVPS = PVPS - BottomEdge.matEvaluateStrongFormEdgeRHS( PVPSBotEFluxM, PVPSBotEFluxP, PVPSBotEFluxS );

PWPSBotEFluxM = BottomEdge.nz .* fm(:,:,11) ./ fm(:,:,4);
PWPSBotEFluxP = BottomEdge.nz .* fp(:,:,11) ./ fp(:,:,4);
PWPSBotEFluxS = ( PWPSBotEFluxM + PWPSBotEFluxP ) ./ 2;
PWPS = PWPS - BottomEdge.matEvaluateStrongFormEdgeRHS( PWPSBotEFluxM, PWPSBotEFluxP, PWPSBotEFluxS );

 [fm, ~] = SurfaceBoundaryEdge.matEvaluateSurfValue( obj.fphys );
 PWPS = PWPS - SurfaceBoundaryEdge.matEvaluateStrongFormEdgeRHS( fm(:,:,11)./fm(:,:,4) .* SurfaceBoundaryEdge.nz, Weta .* SurfaceBoundaryEdge.nz  );


[fm, ~] = BottomBoundaryEdge.matEvaluateSurfValue( obj.fphys );
Wbot = fm(:,:,1)./fm(:,:,4) .* fm(:,:,8) + fm(:,:,2)./fm(:,:,4) .* fm(:,:,9);
PWPS = PWPS - BottomBoundaryEdge.matEvaluateStrongFormEdgeRHS( fm(:,:,11)./fm(:,:,4) .* BottomBoundaryEdge.nz, Wbot .* BottomBoundaryEdge.nz  );

[fm, ~] = BottomBoundaryEdge.matEvaluateSurfValue( obj.fphys );
Unew = fm(:,:,1) ./ fm(:,:,4);
Vnew = fm(:,:,2) ./ fm(:,:,4);
% Wnew = fm(:,:,11) ./ fm(:,:,4);
Wnew = Wbot;

end

function [PUPX, PUPY, PVPX, PVPY, PHPX, PHPY] = CalculateHorizontalDeriveTerm( obj, InnerEdge, BoundaryEdge)
mesh = obj.meshUnion;
cell = mesh.cell;
Hu = obj.fphys{1}(:,:,1);
Hv = obj.fphys{1}(:,:,2);
H = obj.fphys{1}(:,:,4);
u = Hu./H;
v = Hv./H;

PUPX = mesh.rx .* (cell.Dr * u) + mesh.sx .* (cell.Ds * u);
PUPY = mesh.ry .* (cell.Dr * u) + mesh.sy .* (cell.Ds * u);
PVPX = mesh.rx .* (cell.Dr * v) + mesh.sx .* (cell.Ds * v);
PVPY = mesh.ry .* (cell.Dr * v) + mesh.sy .* (cell.Ds * v);
PHPX = mesh.rx .* (cell.Dr * H) + mesh.sx .* (cell.Ds * H);
PHPY = mesh.ry .* (cell.Dr * H) + mesh.sy .* (cell.Ds * H);

[IEfm, IEfp] = InnerEdge.matEvaluateSurfValue( obj.fphys );

PUPXIEfluxM = IEfm(:,:,1)./IEfm(:,:,4) .* InnerEdge.nx;
PUPXIEfluxP = IEfp(:,:,1)./IEfp(:,:,4) .* InnerEdge.nx;
PUPXIEfluxS = ( PUPXIEfluxM +  PUPXIEfluxP ) ./ 2;
PUPX = PUPX - InnerEdge.matEvaluateStrongFormEdgeRHS( PUPXIEfluxM, PUPXIEfluxP, PUPXIEfluxS );

PUPYIEfluxM = IEfm(:,:,1)./IEfm(:,:,4) .* InnerEdge.ny;
PUPYIEfluxP = IEfp(:,:,1)./IEfp(:,:,4) .* InnerEdge.ny;
PUPYIEfluxS = ( PUPYIEfluxM +  PUPYIEfluxP ) ./ 2 ;
PUPY = PUPY - InnerEdge.matEvaluateStrongFormEdgeRHS( PUPYIEfluxM, PUPYIEfluxP, PUPYIEfluxS );

PVPXIEfluxM = IEfm(:,:,2)./IEfm(:,:,4) .* InnerEdge.nx;
PVPXIEfluxP = IEfp(:,:,2)./IEfp(:,:,4) .* InnerEdge.nx;
PVPXIEfluxS = ( PVPXIEfluxM +  PVPXIEfluxP ) ./ 2;
PVPX = PVPX - InnerEdge.matEvaluateStrongFormEdgeRHS( PVPXIEfluxM, PVPXIEfluxP, PVPXIEfluxS );

PVPYIEfluxM = IEfm(:,:,2)./IEfm(:,:,4) .* InnerEdge.ny;
PVPYIEfluxP = IEfp(:,:,2)./IEfp(:,:,4) .* InnerEdge.ny;
PVPYIEfluxS = ( PVPYIEfluxM +  PVPYIEfluxP ) ./ 2 ;
PVPY = PVPY - InnerEdge.matEvaluateStrongFormEdgeRHS( PVPYIEfluxM, PVPYIEfluxP, PVPYIEfluxS );

PHPXIEfluxM = IEfm(:,:,4) .* InnerEdge.nx;
PHPXIEfluxP = IEfp(:,:,4) .* InnerEdge.nx;
PHPXIEfluxS = ( PHPXIEfluxM + PHPXIEfluxP ) ./ 2;
PHPX = PHPX - InnerEdge.matEvaluateStrongFormEdgeRHS( PHPXIEfluxM, PHPXIEfluxP, PHPXIEfluxS );

PHPYIEfluxM = IEfm(:,:,4) .* InnerEdge.ny;
PHPYIEfluxP = IEfp(:,:,4) .* InnerEdge.ny;
PHPYIEfluxS = ( PHPYIEfluxM + PHPYIEfluxP ) ./ 2;
PHPY = PHPY - InnerEdge.matEvaluateStrongFormEdgeRHS( PHPYIEfluxM, PHPYIEfluxP, PHPYIEfluxS );

[BEfm, BEfp] = BoundaryEdge.matEvaluateSurfValue( obj.fphys );
[BEfm, BEfp] = ImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, BEfm, BEfp, obj.fext3d{1} );

PUPXBEfluxM = BEfm(:,:,1)./BEfm(:,:,4) .* BoundaryEdge.nx;
PUPXBEfluxS = ( BEfp(:,:,1)./BEfp(:,:,4) +  BEfm(:,:,1)./BEfm(:,:,4) ) ./ 2 .* BoundaryEdge.nx;
PUPX = PUPX - BoundaryEdge.matEvaluateStrongFormEdgeRHS( PUPXBEfluxM, PUPXBEfluxS );

PUPYBEfluxM = BEfm(:,:,1)./BEfm(:,:,4) .* BoundaryEdge.ny;
PUPYBEfluxS = ( BEfp(:,:,1)./BEfp(:,:,4) +  BEfm(:,:,1)./BEfm(:,:,4) ) ./ 2 .* BoundaryEdge.ny;
PUPY = PUPY - BoundaryEdge.matEvaluateStrongFormEdgeRHS( PUPYBEfluxM, PUPYBEfluxS );

PVPXBEfluxM = BEfm(:,:,2)./BEfm(:,:,4) .* BoundaryEdge.nx;
PVPXBEfluxS = ( BEfp(:,:,2)./BEfp(:,:,4) +  BEfm(:,:,2)./BEfm(:,:,4) ) ./ 2 .* BoundaryEdge.nx;
PVPX = PVPX - BoundaryEdge.matEvaluateStrongFormEdgeRHS( PVPXBEfluxM, PVPXBEfluxS );

PVPYBEfluxM = BEfm(:,:,2)./BEfm(:,:,4) .* BoundaryEdge.ny;
PVPYBEfluxS = ( BEfp(:,:,2)./BEfp(:,:,4) +  BEfm(:,:,2)./BEfm(:,:,4) ) ./ 2 .* BoundaryEdge.ny;
PVPY = PVPY - BoundaryEdge.matEvaluateStrongFormEdgeRHS( PVPYBEfluxM, PVPYBEfluxS );

PHPXBEfluxM = BEfm(:,:,4) .* BoundaryEdge.nx;
PHPXBEfluxS = ( BEfp(:,:,4) + BEfm(:,:,4) ) ./ 2 .* BoundaryEdge.nx;
PHPX = PHPX - BoundaryEdge.matEvaluateStrongFormEdgeRHS( PHPXBEfluxM, PHPXBEfluxS );

PHPYBEfluxM = BEfm(:,:,4) .* BoundaryEdge.ny;
PHPYBEfluxS = ( BEfp(:,:,4) + BEfm(:,:,4) ) ./ 2 .* BoundaryEdge.ny;
PHPY = PHPY - BoundaryEdge.matEvaluateStrongFormEdgeRHS( PHPYBEfluxM, PHPYBEfluxS );
end

function [ fm, fp ] = ImposeBoundaryCondition( edge, nx, ny, fm, fp, fext3d )

ind = ( edge.ftype == enumBoundaryCondition.ClampedDepth );
fp(:, ind, 4) = fext3d(:, ind, 3);

ind = ( edge.ftype == enumBoundaryCondition.ClampedVel );
fp(:, ind, 1) = fext3d(:, ind, 1);
fp(:, ind, 2) = fext3d(:, ind, 2);

ind = ( edge.ftype == enumBoundaryCondition.Clamped );
fp(:, ind, 1) = fext3d(:, ind, 1);
fp(:, ind, 2) = fext3d(:, ind, 2);
fp(:, ind, 4) = fext3d(:, ind, 3);

ind = ( edge.ftype == enumBoundaryCondition.SlipWall );
Hun =  fm( :, ind, 1 ) .* nx(:, ind) + fm( :, ind, 2).* ny(:, ind);
Hvn = -fm( :, ind, 1 ) .* ny(:, ind) + fm( :, ind, 2).* nx(:, ind);
fp(:, ind, 1) = - Hun .* nx(:, ind) - Hvn .* ny(:, ind);
fp(:, ind, 2) = - Hun .* ny(:, ind) + Hvn .* nx(:, ind);
end