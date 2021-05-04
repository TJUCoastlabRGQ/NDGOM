function matPartialDeriveTermTest( obj )

mesh = obj.meshUnion(1);

U = obj.fphys{1}(:,:,1)./obj.fphys{1}(:,:,4);

PUPX =  mesh.rx .* (mesh.cell.Dr *  U ) + ...
    mesh.sx .* ( mesh.cell.Ds *  U );

edge = mesh.InnerEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( obj.fphys );

Um = fm(:,:,1)./fm(:,:,4);
Up = fp(:,:,1)./fp(:,:,4);

% The numerical at inner edge is given as $\hat p = (p^-+p^+)/2$
PUPX = PUPX - edge.matEvaluateStrongFormEdgeRHS( ...
    Um .* edge.nx, Up .* edge.nx , 0.5 * (Um + Up) .* edge.nx);
disp("==========================For PUPX==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(PUPX-obj.NonhydrostaticSolver.PUPX)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(PUPX-obj.NonhydrostaticSolver.PUPX)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.PUPX)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.PUPX)));
disp("=============================================================================");

PUPS =  mesh.rz .* (mesh.cell.Dr *  U ) + ...
    mesh.sz .* ( mesh.cell.Ds *  U ) + ...
    mesh.tz .* ( mesh.cell.Dt *  U );

edge = mesh.BottomEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( obj.fphys );

Um = fm(:,:,1)./fm(:,:,4);
Up = fp(:,:,1)./fp(:,:,4);

% The numerical at inner edge is given as $\hat p = (p^-+p^+)/2$
PUPS = PUPS - edge.matEvaluateStrongFormEdgeRHS( ...
    Um .* edge.nz, Up .* edge.nz , 0.5 * (Um + Up) .* edge.nz);

disp("==========================For PUPS==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(PUPS-obj.NonhydrostaticSolver.PUPS)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(PUPS-obj.NonhydrostaticSolver.PUPS)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.PUPS)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.PUPS)));
disp("=============================================================================");

V = obj.fphys{1}(:,:,2)./obj.fphys{1}(:,:,4);

PVPY =  mesh.ry .* (mesh.cell.Dr *  V ) + ...
    mesh.sy .* ( mesh.cell.Ds *  V );

edge = mesh.InnerEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( obj.fphys );

Vm = fm(:,:,2)./fm(:,:,4);
Vp = fp(:,:,2)./fp(:,:,4);

% The numerical at inner edge is given as $\hat p = (p^-+p^+)/2$
PVPY = PVPY - edge.matEvaluateStrongFormEdgeRHS( ...
    Vm .* edge.ny, Vp .* edge.ny , 0.5 * (Vm + Vp) .* edge.ny);

disp("==========================For PVPY==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(PVPY-obj.NonhydrostaticSolver.PVPY)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(PVPY-obj.NonhydrostaticSolver.PVPY)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.PVPY)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.PVPY)));
disp("=============================================================================");

PVPS =  mesh.rz .* (mesh.cell.Dr *  V ) + ...
    mesh.sz .* ( mesh.cell.Ds *  V ) + ...
    mesh.tz .* ( mesh.cell.Dt *  V );

edge = mesh.BottomEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( obj.fphys );

Vm = fm(:,:,2)./fm(:,:,4);
Vp = fp(:,:,2)./fp(:,:,4);

% The numerical at inner edge is given as $\hat p = (p^-+p^+)/2$
PVPS = PVPS - edge.matEvaluateStrongFormEdgeRHS( ...
    Vm .* edge.nz, Vp .* edge.nz , 0.5 * (Vm + Vp) .* edge.nz);

disp("==========================For PVPS==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(PVPS-obj.NonhydrostaticSolver.PVPS)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(PVPS-obj.NonhydrostaticSolver.PVPS)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.PVPS)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.PVPS)));
disp("=============================================================================");

PSPX = mesh.rx .* (mesh.cell.Dr *  obj.fphys{1}(:,:,4) ) + ...
    mesh.sx .* ( mesh.cell.Ds *  obj.fphys{1}(:,:,4) );

edge = mesh.InnerEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( obj.fphys );
PSPX = PSPX - edge.matEvaluateStrongFormEdgeRHS( ...
    fm(:,:,4) .* edge.nx, fp(:,:,4) .* edge.nx , 0.5 * (fm(:,:,4) + fp(:,:,4)) .* edge.nx);

PSPX = -1 * (1 + mesh.z)./obj.fphys{1}(:,:,4) .* PSPX - 1./obj.fphys{1}(:,:,4) .* obj.fphys{1}(:,:,8);

disp("==========================For PSPX==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(PSPX-obj.NonhydrostaticSolver.PSPX)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(PSPX-obj.NonhydrostaticSolver.PSPX)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.PSPX)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.PSPX)));
disp("=============================================================================");


PSPY = mesh.ry .* (mesh.cell.Dr *  obj.fphys{1}(:,:,4) ) + ...
    mesh.sy .* ( mesh.cell.Ds *  obj.fphys{1}(:,:,4) );

edge = mesh.InnerEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( obj.fphys );
PSPY = PSPY - edge.matEvaluateStrongFormEdgeRHS( ...
    fm(:,:,4) .* edge.ny, fp(:,:,4) .* edge.ny , 0.5 * (fm(:,:,4) + fp(:,:,4)) .* edge.ny);

PSPY = -1 * (1 + mesh.z)./obj.fphys{1}(:,:,4) .* PSPY - 1./obj.fphys{1}(:,:,4) .* obj.fphys{1}(:,:,9);

disp("==========================For PSPY==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(PSPY-obj.NonhydrostaticSolver.PSPY)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(PSPY-obj.NonhydrostaticSolver.PSPY)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.PSPY)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.PSPY)));
disp("=============================================================================");

SQPSPX = PSPX.^2;
SQPSPY = PSPY.^2;

disp("==========================For SQPSPX==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(SQPSPX-obj.NonhydrostaticSolver.SQPSPX)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(SQPSPX-obj.NonhydrostaticSolver.SQPSPX)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.SQPSPX)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.SQPSPX)));
disp("=============================================================================");

disp("==========================For SQPSPY==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(SQPSPY-obj.NonhydrostaticSolver.SQPSPY)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(SQPSPY-obj.NonhydrostaticSolver.SQPSPY)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.SQPSPY)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.SQPSPY)));
disp("=============================================================================");


Data = cell(1);
Data{1} = PSPX;
SPSPX =  mesh.rx .* (mesh.cell.Dr *  PSPX ) + ...
    mesh.sx .* ( mesh.cell.Ds *  PSPX );
[ fm, fp ] = edge.matEvaluateSurfValue( Data );
SPSPX = SPSPX - edge.matEvaluateStrongFormEdgeRHS( ...
    fm .* edge.nx, fp .* edge.nx , 0.5 * (fm + fp) .* edge.nx);

Data{1} = PSPY;
SPSPY =  mesh.ry .* (mesh.cell.Dr *  PSPY ) + ...
    mesh.sy .* ( mesh.cell.Ds *  PSPY );
[ fm, fp ] = edge.matEvaluateSurfValue( Data );
SPSPY = SPSPY - edge.matEvaluateStrongFormEdgeRHS( ...
    fm .* edge.ny, fp .* edge.ny , 0.5 * (fm + fp) .* edge.ny);

disp("==========================For SPSPX==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(SPSPX-obj.NonhydrostaticSolver.SPSPX)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(SPSPX-obj.NonhydrostaticSolver.SPSPX)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.SPSPX)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.SPSPX)));
disp("=============================================================================");

disp("==========================For SPSPY==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(SPSPY-obj.NonhydrostaticSolver.SPSPY)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(SPSPY-obj.NonhydrostaticSolver.SPSPY)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.SPSPY)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.SPSPY)));
disp("=============================================================================");









W = cell(1);
W{1} = obj.fphys{1}(:,:,11)./obj.fphys{1}(:,:,4);

%> Volume integral
PWPS =  mesh.rz .* (mesh.cell.Dr *  W{1} ) + ...
    mesh.sz .* ( mesh.cell.Ds *  W{1} ) + ...
    mesh.tz .* ( mesh.cell.Dt *  W{1} );

Etax = mesh.rx .* (mesh.cell.Dr * obj.fphys{1}(:,:,4)) + ...
    mesh.sx .* ( mesh.cell.Ds * obj.fphys{1}(:,:,4) ) + ...
    mesh.tx .* ( mesh.cell.Dt * obj.fphys{1}(:,:,4) );

Etay = mesh.ry .* (mesh.cell.Dr * obj.fphys{1}(:,:,4)) + ...
    mesh.sy .* ( mesh.cell.Ds * obj.fphys{1}(:,:,4) ) + ...
    mesh.ty .* ( mesh.cell.Dt * obj.fphys{1}(:,:,4) );

edge = mesh.InnerEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( obj.fphys );

% The numerical at inner edge is given as $\hat p = (p^-+p^+)/2$
Etax = Etax - edge.matEvaluateStrongFormEdgeRHS( ...
    fm(:,:,4) .* edge.nx, fp(:,:,4) .* edge.nx , 0.5 * (fm(:,:,4) + fp(:,:,4)) .* edge.nx);

Etay = Etay - edge.matEvaluateStrongFormEdgeRHS( ...
    fm(:,:,4) .* edge.ny, fp(:,:,4) .* edge.ny , 0.5 * (fm(:,:,4) + fp(:,:,4)) .* edge.ny);

edge = mesh.BoundaryEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( obj.fphys );

Etax = Etax - edge.matEvaluateStrongFormEdgeRHS( ...
    fm(:,:,4) .* edge.nx, 0.5 * (fm(:,:,4) + fp(:,:,4)) .* edge.nx) + obj.fphys{1}(:,:,8);

Etay = Etay - edge.matEvaluateStrongFormEdgeRHS( ...
    fm(:,:,4) .* edge.ny, 0.5 * (fm(:,:,4) + fp(:,:,4)) .* edge.ny) + obj.fphys{1}(:,:,9);

edge = mesh.BottomEdge;
[ fm, fp ] = edge.matEvaluateSurfValue( W );

PWPS = PWPS - edge.matEvaluateStrongFormEdgeRHS( ...
    fm .* edge.nz, fp .* edge.nz , 0.5 * (fm + fp) .* edge.nz);

edge = mesh.BottomBoundaryEdge;
[fm, ~] = edge.matEvaluateSurfValue( obj.fphys );
Wbot = fm(:,:,1)./fm(:,:,4).*fm(:,:,8) + fm(:,:,2)./fm(:,:,4).*fm(:,:,9);
PWPS = PWPS - edge.matEvaluateStrongFormEdgeRHS( ...
    fm(:,:,11)./fm(:,:,4) .* edge.nz, Wbot .* edge.nz);

edge = mesh.SurfaceBoundaryEdge;
[fm, ~] = edge.matEvaluateSurfValue( obj.fphys );
Usurf = fm(:,:,1)./fm(:,:,4);
Vsurf = fm(:,:,2)./fm(:,:,4);


Etat2d = -1 * ( mesh.mesh2d.rx .* (mesh.mesh2d.cell.Dr * obj.fphys2d{1}(:,:,2)) + ...
    mesh.mesh2d.sx .* (mesh.mesh2d.cell.Ds * obj.fphys2d{1}(:,:,2)) + ...
    mesh.mesh2d.ry .* (mesh.mesh2d.cell.Dr * obj.fphys2d{1}(:,:,3)) + ...
    mesh.mesh2d.sy .* (mesh.mesh2d.cell.Ds * obj.fphys2d{1}(:,:,3)) );

edge = mesh.mesh2d.InnerEdge;
[fm, fp] = edge.matEvaluateSurfValue( obj.fphys2d );
Inputfm = zeros(edge.Nfp, edge.Ne, 3);
Inputfp = zeros(edge.Nfp, edge.Ne, 3);
Inputfm(:,:,1) = fm(:,:,2);Inputfm(:,:,2) = fm(:,:,3);Inputfm(:,:,3) = fm(:,:,1);
Inputfp(:,:,1) = fp(:,:,2);Inputfp(:,:,2) = fp(:,:,3);Inputfp(:,:,3) = fp(:,:,1);
FluxS = mxCalculateInnerEdgeNumFluxTerm( struct(edge), Inputfm, Inputfp, obj.gra, obj.hcrit);

Etat2d = Etat2d + edge.matEvaluateStrongFormEdgeRHS( fm(:,:,2).*edge.nx + fm(:,:,3).*edge.ny,...
    fp(:,:,2).*edge.nx + fp(:,:,3).*edge.ny, FluxS);

edge = mesh.mesh2d.BoundaryEdge;
[fm, fp] = edge.matEvaluateSurfValue( obj.fphys2d );
Inputfm = zeros(edge.Nfp, edge.Ne, 3);
Inputfp = zeros(edge.Nfp, edge.Ne, 3);
Inputfm(:,:,1) = fm(:,:,2);Inputfm(:,:,2) = fm(:,:,3);Inputfm(:,:,3) = fm(:,:,1);
Inputfp(:,:,1) = fp(:,:,2);Inputfp(:,:,2) = fp(:,:,3);Inputfp(:,:,3) = fp(:,:,1);
zM = fm(:,:,4); zP = fp(:,:,4);

FluxS = mxCalculateBoundaryEdgeNumFluxTerm( struct(edge), Inputfm, Inputfp, obj.gra, obj.hcrit, zM, zP, edge.ftype, obj.fext2d{1});

Etat2d = Etat2d + edge.matEvaluateStrongFormEdgeRHS( fm(:,:,2).*edge.nx + fm(:,:,3).*edge.ny, FluxS);


TempEtax = cell(1); TempEtax{1} = Etax;
TempEtay = cell(1); TempEtay{1} = Etay;

edge = mesh.SurfaceBoundaryEdge;
[EtaSurfx, ~] = edge.matEvaluateSurfValue( TempEtax );
[EtaSurfy, ~] = edge.matEvaluateSurfValue( TempEtay );

Wsurf = Etat2d + Usurf .* EtaSurfx + Vsurf .* EtaSurfy;

[fm, ~] = edge.matEvaluateSurfValue( obj.fphys );

 PWPS = PWPS - edge.matEvaluateStrongFormEdgeRHS( ...
     fm(:,:,11)./fm(:,:,4) .* edge.nz, Wsurf .* edge.nz);

disp("==========================For PWPS==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(PWPS-obj.NonhydrostaticSolver.PWPS)));
disp("The minimum difference between the matlab version and C version is:");
disp(min(min(PWPS-obj.NonhydrostaticSolver.PWPS)));
disp("The maximum value of the C version is:");
disp(max(max(obj.NonhydrostaticSolver.PWPS)));
disp("The minimum value of the C version is:");
disp(min(min(obj.NonhydrostaticSolver.PWPS)));
disp("=============================================================================");

end