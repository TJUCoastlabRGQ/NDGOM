function matfphys = matUpdateConservativeFinalVelocity( obj, matfphys, NonhydroPressure, dt, PSPX, PSPY )
%MATUPDATECONSERVATIVEFINALVELOCITY 此处显示有关此函数的摘要
%   此处显示详细说明
mesh = obj.meshUnion(1);
cell = mesh.cell;
qx = mesh.rx .* (cell.Dr * NonhydroPressure{1}) + mesh.sx .* (cell.Ds * NonhydroPressure{1});
qy = mesh.ry .* (cell.Dr * NonhydroPressure{1}) + mesh.sy .* (cell.Ds * NonhydroPressure{1});
qz = mesh.tz .* (cell.Dt * NonhydroPressure{1});
edge = obj.meshUnion.InnerEdge;
[Npm, Npp] = edge.matEvaluateSurfValue(NonhydroPressure);

pxIEfluxM = Npm .* edge.nx;
pxIEfluxP = Npp .* edge.nx;
pxIEfluxS = (pxIEfluxM + pxIEfluxP)./2;
qx = qx - edge.matEvaluateStrongFormEdgeRHS( pxIEfluxM, pxIEfluxP, pxIEfluxS );

pyIEfluxM = Npm .* edge.ny;
pyIEfluxP = Npp .* edge.ny;
pyIEfluxS = (pyIEfluxM + pyIEfluxP)./2;
qy = qy - edge.matEvaluateStrongFormEdgeRHS( pyIEfluxM, pyIEfluxP, pyIEfluxS );

edge = obj.meshUnion.BottomEdge; 
[Npm, Npp] = edge.matEvaluateSurfValue(NonhydroPressure);
pzBotEfluxM = Npm .* edge.nz;
pzBotEfluxP = Npp .* edge.nz;
pzBotEfluxS = (pzBotEfluxM + pzBotEfluxP)./2;
qz = qz - edge.matEvaluateStrongFormEdgeRHS( pzBotEfluxM, pzBotEfluxP, pzBotEfluxS );

edge = obj.meshUnion.SurfaceBoundaryEdge; 
[Npm, ~] = edge.matEvaluateSurfValue(NonhydroPressure);
pzSurfBEfluxM = Npm .* edge.nz;
pzSurfBEfluxS = 0 .* pzSurfBEfluxM;
qz = qz - edge.matEvaluateStrongFormEdgeRHS( pzSurfBEfluxM, pzSurfBEfluxS );

matfphys(:,:,1) = matfphys(:,:,1) - dt*matfphys(:,:,4)/obj.NonhydrostaticSolver.rho.*(qx + qz .* PSPX);
matfphys(:,:,2) = matfphys(:,:,2) - dt*matfphys(:,:,4)/obj.NonhydrostaticSolver.rho.*(qy + qz .* PSPY);
matfphys(:,:,11) = matfphys(:,:,11) - dt/obj.NonhydrostaticSolver.rho * qz;
end

