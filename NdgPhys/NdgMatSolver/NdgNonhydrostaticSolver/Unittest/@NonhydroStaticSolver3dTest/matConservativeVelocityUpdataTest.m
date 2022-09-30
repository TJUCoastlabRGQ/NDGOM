function matConservativeVelocityUpdataTest( obj )
Nonhydro = cell(1);
Nonhydro{1} = 1000 * rand( obj.meshUnion.cell.Np, obj.meshUnion.K );

mesh = obj.meshUnion(1);

PNPX =  mesh.rx .* (mesh.cell.Dr *  Nonhydro{1} ) + ...
    mesh.sx .* ( mesh.cell.Ds *  Nonhydro{1} );

edge = mesh.InnerEdge;
[ Pm, Pp ] = edge.matEvaluateSurfValue( Nonhydro );

PNPX = PNPX - edge.matEvaluateStrongFormEdgeRHS( ...
    Pm .* edge.nx, Pp .* edge.nx , 0.5 * (Pm + Pp) .* edge.nx);

PNPY =  mesh.ry .* (mesh.cell.Dr *  Nonhydro{1} ) + ...
    mesh.sy .* ( mesh.cell.Ds *  Nonhydro{1} );

edge = mesh.InnerEdge;
[ Pm, Pp ] = edge.matEvaluateSurfValue( Nonhydro );

PNPY = PNPY - edge.matEvaluateStrongFormEdgeRHS( ...
    Pm .* edge.ny, Pp .* edge.ny , 0.5 * (Pm + Pp) .* edge.ny);

PNPS =  mesh.tz .* (mesh.cell.Dt *  Nonhydro{1} );

edge = mesh.BottomEdge;
[ Pm, Pp ] = edge.matEvaluateSurfValue( Nonhydro );

PNPS = PNPS - edge.matEvaluateStrongFormEdgeRHS( ...
    Pm .* edge.nz, Pp .* edge.nz , 0.5 * (Pm + Pp) .* edge.nz);

edge = mesh.SurfaceBoundaryEdge;
[ Pm, ~ ] = edge.matEvaluateSurfValue( Nonhydro );

PNPS = PNPS - edge.matEvaluateStrongFormEdgeRHS( ...
    Pm .* edge.nz, zeros(size(edge.nz)));

dt = 0.052;

Hu = obj.fphys{1}(:,:,1) - dt * obj.fphys{1}(:,:,4)/1000 .* ( PNPX + PNPS.*obj.NonhydrostaticSolver.PSPX );
Hv = obj.fphys{1}(:,:,2) - dt * obj.fphys{1}(:,:,4)/1000 .* ( PNPY + PNPS.*obj.NonhydrostaticSolver.PSPY );
Hw = obj.fphys{1}(:,:,11) - dt/1000 .* ( PNPS );

Data = obj.NonhydrostaticSolver.TestUpdateConservativeFinalVelocity( obj, Nonhydro{1}, obj.fphys, obj.fphys2d, dt );

disp("==========================For Hu==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(Hu-Data{1}(:,:,1))));
disp("=============================================================================");

disp("==========================For Hv==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(Hv-Data{1}(:,:,2))));
disp("=============================================================================");

disp("==========================For Hw==========================================");
disp("The maximum difference between the matlab version and C version is:");
disp(max(max(Hw-Data{1}(:,:,11))));
disp("=============================================================================");

end