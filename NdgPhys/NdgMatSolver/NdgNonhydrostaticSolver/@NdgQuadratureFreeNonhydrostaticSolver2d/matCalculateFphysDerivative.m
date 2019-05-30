function matCalculateFphysDerivative(obj, mesh, fphys, physClass)
InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
% Inner edge part for $\frac{(H + 2*b)}{x}$
[ fluxMX_HB, fluxPX_HB ] = getLocalAndAdjacentFluxTerm( InnerEdge, InnerEdge.nx, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]));
fluxX_HB = getCentralFluxTerm( InnerEdge, InnerEdge.nx, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]) );
fluxX_HB = getUpwindFluxTerm( InnerEdge, InnerEdge.nx, fm, fp, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), fluxX_HB );
obj.HBx = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX_HB, fluxPX_HB, fluxX_HB);
% Inner edge part for $\frac{(H + 2*b)}{y}$
[ fluxMY_HB, fluxPY_HB ] = getLocalAndAdjacentFluxTerm( InnerEdge, InnerEdge.ny, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]));
fluxY_HB = getCentralFluxTerm( InnerEdge, InnerEdge.ny, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]) );
fluxY_HB = getUpwindFluxTerm( InnerEdge, InnerEdge.ny, fm, fp, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), fluxY_HB );
obj.HBy = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY_HB, fluxPY_HB, fluxY_HB);
% Inner edge part for $\frac{H}{x}$  
[ fluxMX_H, fluxPX_H ] = getLocalAndAdjacentFluxTerm( InnerEdge, InnerEdge.nx, num2cell(fphys{1}(:,:,1), [1 2]));
fluxX_H = getCentralFluxTerm( InnerEdge, InnerEdge.nx, num2cell(fphys{1}(:,:,1), [1 2]) );
fluxX_H = getUpwindFluxTerm( InnerEdge, InnerEdge.nx, fm, fp, num2cell(fphys{1}(:,:,1),[1 2]), fluxX_H );
obj.fhx = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX_H, fluxPX_H, fluxX_H);
% Inner edge part for $\frac{H}{y}$  
[ fluxMY_H, fluxPY_H ] = getLocalAndAdjacentFluxTerm( InnerEdge, InnerEdge.ny, num2cell(fphys{1}(:,:,1), [1 2]));
fluxY_H = getCentralFluxTerm( InnerEdge, InnerEdge.ny, num2cell(fphys{1}(:,:,1), [1 2]) );
fluxY_H = getUpwindFluxTerm( InnerEdge, InnerEdge.ny, fm, fp, num2cell(fphys{1}(:,:,1),[1 2]), fluxY_H );
obj.fhy = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY_H, fluxPY_H, fluxY_H);
% Inner edge part for $\frac{Hu}{x}$  
[ fluxMX_Hu, fluxPX_Hu ] = getLocalAndAdjacentFluxTerm( InnerEdge, InnerEdge.nx, num2cell(fphys{1}(:,:,2), [1 2]));
fluxX_Hu = getCentralFluxTerm( InnerEdge, InnerEdge.nx, num2cell(fphys{1}(:,:,2), [1 2]) );
fluxX_Hu = getUpwindFluxTerm( InnerEdge, InnerEdge.nx, fm, fp, num2cell(fphys{1}(:,:,2),[1 2]), fluxX_Hu );
obj.hux = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX_Hu, fluxPX_Hu, fluxX_Hu);
% Inner edge part for $\frac{Hv}{y}$  
[ fluxMY_Hv, fluxPY_Hv ] = getLocalAndAdjacentFluxTerm( InnerEdge, InnerEdge.ny, num2cell(fphys{1}(:,:,3), [1 2]));
fluxY_Hv = getCentralFluxTerm( InnerEdge, InnerEdge.ny, num2cell(fphys{1}(:,:,3), [1 2]) );
fluxY_Hv = getUpwindFluxTerm( InnerEdge, InnerEdge.ny, fm, fp, num2cell(fphys{1}(:,:,3),[1 2]), fluxY_Hv );
obj.hvy = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY_Hv, fluxPY_Hv, fluxY_Hv);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys );
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1} );

% Boundary edge part for $\frac{(H + 2*b)}{x}$
[ fluxMX_HB, ~ ] = getLocalAndAdjacentFluxTerm( BoundaryEdge, BoundaryEdge.nx, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]));
fluxX_HB = getCentralFluxTerm( BoundaryEdge, BoundaryEdge.nx, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]) );
fluxX_HB = getUpwindFluxTerm( BoundaryEdge, BoundaryEdge.nx, fm, fp, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), fluxX_HB );
obj.HBx = -obj.HBx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX_HB, fluxX_HB);
% Boundary edge part for $\frac{(H + 2*b)}{y}$
[ fluxMY_HB, ~ ] = getLocalAndAdjacentFluxTerm( BoundaryEdge, BoundaryEdge.ny, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]));
fluxY_HB = getCentralFluxTerm( BoundaryEdge, BoundaryEdge.ny, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]) );
fluxY_HB = getUpwindFluxTerm( BoundaryEdge, BoundaryEdge.ny, fm, fp, num2cell(fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4),[1 2]), fluxY_HB );
obj.HBy = -obj.HBy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY_HB, fluxY_HB);
% Boundary edge part for $\frac{H}{x}$  
[ fluxMX_H, ~ ] = getLocalAndAdjacentFluxTerm( BoundaryEdge, BoundaryEdge.nx, num2cell(fphys{1}(:,:,1), [1 2]));
fluxX_H = getCentralFluxTerm( BoundaryEdge, BoundaryEdge.nx, num2cell(fphys{1}(:,:,1), [1 2]) );
fluxX_H = getUpwindFluxTerm( BoundaryEdge, BoundaryEdge.nx, fm, fp, num2cell(fphys{1}(:,:,1),[1 2]), fluxX_H );
obj.fhx = -obj.fhx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX_H, fluxX_H);
% Boundary edge part for $\frac{H}{y}$  
[ fluxMY_H, ~ ] = getLocalAndAdjacentFluxTerm( BoundaryEdge, BoundaryEdge.ny, num2cell(fphys{1}(:,:,1), [1 2]));
fluxY_H = getCentralFluxTerm( BoundaryEdge, BoundaryEdge.ny, num2cell(fphys{1}(:,:,1), [1 2]) );
fluxY_H = getUpwindFluxTerm( BoundaryEdge, BoundaryEdge.ny, fm, fp, num2cell(fphys{1}(:,:,1),[1 2]), fluxY_H );
obj.fhy = -obj.fhy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY_H, fluxY_H);
% Boundary edge part for $\frac{Hu}{x}$  
[ fluxMX_Hu, ~ ] = getLocalAndAdjacentFluxTerm( BoundaryEdge, BoundaryEdge.nx, num2cell(fphys{1}(:,:,2), [1 2]));
fluxX_Hu = getCentralFluxTerm( BoundaryEdge, BoundaryEdge.nx, num2cell(fphys{1}(:,:,2), [1 2]) );
fluxX_Hu = getUpwindFluxTerm( BoundaryEdge, BoundaryEdge.nx, fm, fp, num2cell(fphys{1}(:,:,2),[1 2]), fluxX_Hu );
obj.hux = -obj.hux - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX_Hu, zeros(size(fluxX_Hu)));
% Boundary edge part for $\frac{Hv}{y}$  
[ fluxMY_Hv, ~ ] = getLocalAndAdjacentFluxTerm( BoundaryEdge, BoundaryEdge.ny, num2cell(fphys{1}(:,:,3), [1 2]));
fluxY_Hv = getCentralFluxTerm( BoundaryEdge, BoundaryEdge.ny, num2cell(fphys{1}(:,:,3), [1 2]) );
fluxY_Hv = getUpwindFluxTerm( BoundaryEdge, BoundaryEdge.ny, fm, fp, num2cell(fphys{1}(:,:,3),[1 2]), fluxY_Hv );
obj.hvy = -obj.hvy - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY_Hv, zeros(size(fluxY_Hv)));

% Volume part for $\frac{(H + 2*b)}{x}$ and $\frac{(H + 2*b)}{y}$ 
[VHBX, VHBY] = obj.matVolumeIntegral( mesh, fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4), fphys{1}(:,:,1) + 2 * fphys{1}(:,:,4));
obj.HBx = obj.HBx + VHBX;
obj.HBy = obj.HBy + VHBY;
% Volume part for $\frac{H}{x}$ and $\frac{H}{y}$ 
[VHX, VHY] = obj.matVolumeIntegral( mesh, fphys{1}(:,:,1), fphys{1}(:,:,1));
obj.fhx = obj.fhx + VHX;
obj.fhy = obj.fhy + VHY;
% Volume part for $\frac{Hu}{x}$
[VHUX, VHVY] = obj.matVolumeIntegral( mesh, fphys{1}(:,:,2), fphys{1}(:,:,3));
obj.hux = obj.hux + VHUX;
% Volume part for $\frac{Hv}{y}$
obj.hvy = obj.hvy + VHVY;

% Calculation for H2Bx and H2By
% Inner edge part for $\frac{\frac{H+2b}{x}}{x}$ 
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[ fluxMX_H2Bx, fluxPX_H2Bx ] = getLocalAndAdjacentFluxTerm( InnerEdge, InnerEdge.nx, num2cell(obj.HBx,[1 2]));
fluxX_H2Bx = getCentralFluxTerm( InnerEdge, InnerEdge.nx, num2cell(obj.HBx,[1 2]) );
fluxX_H2Bx = getDownwindFluxTerm( InnerEdge, InnerEdge.nx, fm, fp, num2cell(obj.HBx,[1 2]), fluxX_H2Bx );
obj.H2Bx = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMX_H2Bx, fluxPX_H2Bx, fluxX_H2Bx);
% Inner edge part for $\frac{\frac{H+2b}{y}}{y}$ 
[ fluxMY_H2By, fluxPY_H2By ] = getLocalAndAdjacentFluxTerm( InnerEdge, InnerEdge.ny, num2cell(obj.HBy,[1 2]));
fluxY_H2By = getCentralFluxTerm( InnerEdge, InnerEdge.ny, num2cell(obj.HBy,[1 2]) );
fluxY_H2By = getDownwindFluxTerm( InnerEdge, InnerEdge.ny, fm, fp, num2cell(obj.HBy,[1 2]), fluxY_H2By );
obj.H2By = InnerEdge.matEvaluateStrongFromEdgeRHS(fluxMY_H2By, fluxPY_H2By, fluxY_H2By);

% Boundary edge part for $\frac{\frac{H+2b}{x}}{x}$  
[ fluxMX_HBx, ~ ] = getLocalAndAdjacentFluxTerm( BoundaryEdge, BoundaryEdge.nx, num2cell(obj.HBx,[1 2]));
obj.H2Bx = -obj.H2Bx - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMX_HBx, zeros(size(fluxMX_HBx)) );
% Boundary edge part for $\frac{\frac{H+2b}{y}}{y}$  
[ fluxMY_HBy, ~ ] = getLocalAndAdjacentFluxTerm( BoundaryEdge, BoundaryEdge.ny, num2cell(obj.HBy,[1 2]));
obj.H2By = -obj.H2By - BoundaryEdge.matEvaluateStrongFromEdgeRHS(fluxMY_HBy, zeros(size(fluxMY_HBy)) );

% Volume part for $\frac{\frac{H+2b}{x}}{x}$  
[VH2BxX, VH2ByY] = obj.matVolumeIntegral( mesh, obj.HBx, obj.HBy);
obj.H2Bx = obj.H2Bx + VH2BxX;
% Volume part for $\frac{\frac{H+2b}{y}}{y}$
obj.H2By = obj.H2By + VH2ByY;
end

function [fluxM, fluxP] = getLocalAndAdjacentFluxTerm(edge, vector, Variable)
[vfm, vfp] = edge.matEvaluateSurfValue(Variable);
fluxM = vector .* vfm;
fluxP = vector .* vfp;
end

function flux = getCentralFluxTerm(edge, vector, Variable)
[vfm, vfp] = edge.matEvaluateSurfValue(Variable);
flux = vector .* ( vfm + vfp )./2;
end

function flux = getUpwindFluxTerm(edge, vector, fm, fp, Variable, flux)
temphum = fm(:,:,2); temphvm = fm(:,:,3);
temphup = fp(:,:,2); temphvp = fp(:,:,3);
[vfm, vfp] = edge.matEvaluateSurfValue(Variable);

Index = ( temphum .* edge.nx + temphvm .* edge.ny > 0 & - temphup .* edge.nx - temphvp .* edge.ny <= 0 );
flux(Index) = vfm(Index) .* vector(Index);

Index = ( temphum .* edge.nx + temphvm .* edge.ny <= 0 & - temphup .* edge.nx - temphvp .* edge.ny > 0 );
flux( Index ) =  vfp(Index) .* vector(Index);
end

function [flux] = getDownwindFluxTerm(edge, vector, fm, fp, Variable, flux)
temphum = fm(:,:,2); temphvm = fm(:,:,3);
temphup = fp(:,:,2); temphvp = fp(:,:,3);
[vfm, vfp] = edge.matEvaluateSurfValue(Variable);

Index = ( temphum .* edge.nx + temphvm .* edge.ny > 0 & - temphup .* edge.nx - temphvp .* edge.ny <= 0 );
flux(Index) = vfp(Index) .* vector(Index);

Index = ( temphum .* edge.nx + temphvm .* edge.ny <= 0 & - temphup .* edge.nx - temphvp .* edge.ny > 0 );
flux( Index ) =  vfm(Index) .* vector(Index);
end