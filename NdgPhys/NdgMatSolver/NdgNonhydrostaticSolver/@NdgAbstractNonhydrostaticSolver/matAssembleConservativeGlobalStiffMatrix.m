function [ StiffMatrix ] = matAssembleConservativeGlobalStiffMatrix( obj, hm, hp, height, physClass )
mesh = physClass.meshUnion(1);
K = mesh.K; Np = mesh.cell.Np;
StiffMatrix = zeros( K*Np ); gmat = zeros( Np, K );
% gmat1 = zeros( Np, K );

Nonhydro = zeros( Np, K );
rho = physClass.rho;
dt = physClass.dt;
bx = physClass.zGrad{1}(:,:,1);
by = physClass.zGrad{1}(:,:,2);
for i = 1:K*Np
    gmat(i) = 1.0/2.0 * height(i);
%     gmat1(i) = 1;
    
    Nonhydro(i) = 1;
    fm = gmat(mesh.eidM); fp =   1.0/2.0*hp.*Nonhydro(mesh.eidP);
%     fm1 = gmat1(mesh.eidM);  fp1 = - gmat1(mesh.eidP);
    
    varflux = obj.evaluateVarSurfNumFlux(fm, fp);
%     varflux1 = obj.evaluateVarSurfNumFlux(fm1,fp1);
    
    deltaVarFlux = ( -fm ) - ( -varflux );
%     deltaVarFlux1 = (-fm1) - ( -varflux1 );
    
    [qx, qy] = obj.evaluateAuxialaryVariable(1, gmat, deltaVarFlux );
%     [qx1, qy1] = obj.evaluateAuxialaryVariable(1,gmat1, deltaVarFlux1);
%     df = fm -fp;
    
    [qxm, qxp] = obj.matEvaluateNonhydrostaticSurfaceValue(mesh, qx);
%     [qxm1, qxp1] = obj.matEvaluateNonhydrostaticSurfaceValue(mesh, qx1);   
    
    [qym, qyp] = obj.matEvaluateNonhydrostaticSurfaceValue(mesh, qy);
%     [qym1, qyp1] = obj.matEvaluateNonhydrostaticSurfaceValue(mesh, qy1);
    
    penaltyqx = zeros(size(mesh.eidM));penaltyqy = zeros(size(mesh.eidM));
%     penaltyqx1 = zeros(size(mesh.eidM));penaltyqy1 = zeros(size(mesh.eidM));
    
    fluxqx = obj.matEvaluateNonconservativeFlux(  qxm, qxp, penaltyqx);
%     fluxqx1 = obj.matEvaluateNonconservativeFlux(  qxm1, qxp1, penaltyqx1);
    
    fluxqy = obj.matEvaluateNonconservativeFlux(  qym, qyp, penaltyqy);
%     fluxqy1 = obj.matEvaluateNonconservativeFlux(  qym1, qyp1, penaltyqy1);
    
    divqx = obj.matEvaluateDivergenceDirectionX(mesh, qx);
%     divqx1 = obj.matEvaluateDivergenceDirectionX(mesh, qx1);
    
    divqy = obj.matEvaluateDivergenceDirectionY(mesh, qy);
%     divqy1 = obj.matEvaluateDivergenceDirectionY(mesh, qy1);
    
    deltafluxqx = obj.matEvaluateDeltaSurfaceFluxX( mesh, qxm, fluxqx );
%     deltafluxqx1 = obj.matEvaluateDeltaSurfaceFluxX( mesh, qxm1, fluxqx1 ); 
    
    deltafluxqy = obj.matEvaluateDeltaSurfaceFluxY( mesh, qym, fluxqy );
%     deltafluxqy1 = obj.matEvaluateDeltaSurfaceFluxY( mesh, qym1, fluxqy1 );
    
%     DATAX =  mesh.J .* (mesh.cell.M *((divqx - mesh.cell.LIFT*(mesh.Js .* deltafluxqx )./mesh.J)));
%     DATAX1 =  mesh.J .* (mesh.cell.M *((divqx1 - mesh.cell.LIFT*(mesh.Js .* deltafluxqx1 )./mesh.J)));
    
%     DATAY =  mesh.J .* (mesh.cell.M *((divqy - mesh.cell.LIFT*(mesh.Js .* deltafluxqy )./mesh.J))); 
%     DATAY1 =  mesh.J .* (mesh.cell.M *((divqy1 - mesh.cell.LIFT*(mesh.Js .* deltafluxqy1 )./mesh.J))); 
    
    qbxm = Nonhydro(mesh.eidM) .* bx(mesh.eidM);qbxp = Nonhydro(mesh.eidP) .* bx(mesh.eidP);
    qbym = Nonhydro(mesh.eidM) .* by(mesh.eidM);qbyp = Nonhydro(mesh.eidP) .* by(mesh.eidP);
    penaltyqbx  =  zeros(size(qbxm));penaltyqby  =  zeros(size(qbym));
    fluxqbx = obj.matEvaluateNonconservativeFlux(  qbxm, qbxp, penaltyqbx);
    fluxqby = obj.matEvaluateNonconservativeFlux(  qbym, qbyp, penaltyqby);    
    divqbx = obj.matEvaluateDivergenceDirectionX( mesh, Nonhydro.*bx );
    divqby = obj.matEvaluateDivergenceDirectionY( mesh, Nonhydro.*by );
    deltafluxqbx = obj.matEvaluateDeltaSurfaceFluxX( mesh, Nonhydro( mesh.eidM ) .* bx( mesh.eidM ), fluxqbx );
    deltafluxqby = obj.matEvaluateDeltaSurfaceFluxY( mesh, Nonhydro( mesh.eidM ) .* by( mesh.eidM ), fluxqby );
    
    penaltyhx = zeros(size(mesh.eidM)); penaltyhy = zeros(size(mesh.eidM));
    fluxhx = obj.matEvaluateNonconservativeFlux( hm, hp, penaltyhx);   
    fluxhy = obj.matEvaluateNonconservativeFlux( hm, hp, penaltyhy);
    divhx = obj.matEvaluateDivergenceDirectionX( mesh, height);
    divhy = obj.matEvaluateDivergenceDirectionY( mesh, height);
    deltafluxhx = obj.matEvaluateDeltaSurfaceFluxX( mesh, hm, fluxhx );
    deltafluxhy = obj.matEvaluateDeltaSurfaceFluxY( mesh, hm, fluxhy ); 
    
    rhsu = 2 * dt/rho * mesh.J .* (mesh.cell.M * Nonhydro) + ...
        2 * dt/rho * mesh.J .* (mesh.cell.M * ((qx + Nonhydro .* bx ) .* bx)) + ...
        2 * dt/rho * mesh.J .* (mesh.cell.M * ((qy + Nonhydro .* by ) .* by)) -  ...
        dt/rho * mesh.J .* (mesh.cell.M *((divqx - mesh.cell.LIFT*(mesh.Js .* deltafluxqx )./mesh.J).* height)) - ...
        dt/rho * mesh.J .* (mesh.cell.M *((divqy - mesh.cell.LIFT*(mesh.Js .* deltafluxqy )./mesh.J).* height)) - ...
        dt/rho * mesh.J .* (mesh.cell.M *((divqbx - mesh.cell.LIFT*(mesh.Js .* deltafluxqbx)./mesh.J).* height)) - ...
        dt/rho * mesh.J .* (mesh.cell.M *((divqby - mesh.cell.LIFT*(mesh.Js .* deltafluxqby)./mesh.J).* height))  + ...
        dt/rho * mesh.J .* (mesh.cell.M *((divhx - mesh.cell.LIFT*(mesh.Js .* deltafluxhx)./mesh.J).* qx)) + ...
        dt/rho * mesh.J .* (mesh.cell.M *((divhy - mesh.cell.LIFT*(mesh.Js .* deltafluxhy)./mesh.J).* qy)) + ...
        dt/rho * mesh.J .* (mesh.cell.M *((divhx - mesh.cell.LIFT*(mesh.Js .* deltafluxhx)./mesh.J).* ( Nonhydro .* bx ))) + ...
        dt/rho * mesh.J .* (mesh.cell.M *((divhy - mesh.cell.LIFT*(mesh.Js .* deltafluxhy)./mesh.J).* ( Nonhydro .* by ))) ;
    
    
    StiffMatrix(:,i) = rhsu(:);
    gmat(i) = 0;
    gmat1(i) = 0;
    Nonhydro(i) = 0;
end
end

