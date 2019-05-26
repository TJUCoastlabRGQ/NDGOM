function fphys = matUpdateConservativeFinalVelocity(obj, NonhydroPre , physClass, fphys)
mesh = physClass.meshUnion(1);
% InnerEdge = mesh.InnerEdge;
% BoundaryEdge = mesh.BoundaryEdge;
NonhydrostaticPressure = zeros(mesh.cell.Np, mesh.K);
NonhydroPre = reshape(NonhydroPre, mesh.cell.Np, obj.WetNum);
for i = 1:numel(obj.WetCellIndex)
    NonhydrostaticPressure(:,obj.WetCellIndex(i)) = NonhydroPre(:,i);
end

% NonhydroVolumeflux = 1/2 * NonhydrostaticPressure .* fphys{1}(:,:,1);

% [ NqHx , NqHy ] = obj.matCalculateCharacteristicMatrix( mesh,  BoundaryEdge, InnerEdge, num2cell(NonhydroVolumeflux,[1 2]), num2cell(NonhydroVolumeflux,[1 2]), enumNonhydroBoundaryCondition.Zero);
[ NqHx , NqHy ] = getNonhydroPartialCorrectionTerm( obj, mesh, fphys{1}(:,:,1), NonhydroPre );

% fphys{1}(:,:,6) = fphys{1}(:,:,6) + obj.dt/physClass.rho .* NonhydrostaticPressure;
% fphys{1}(:,:,2) = fphys{1}(:,:,2) - obj.dt/physClass.rho .* (NqHx +NonhydrostaticPressure.*obj.bx);
% fphys{1}(:,:,3) = fphys{1}(:,:,3) - obj.dt/physClass.rho .* (NqHy +NonhydrostaticPressure.*obj.by);

fphys{1}(:,:,6) = fphys{1}(:,:,6) + obj.dt .* NonhydrostaticPressure;
fphys{1}(:,:,2) = fphys{1}(:,:,2) - 1/2 * obj.dt .* ( NqHx );
fphys{1}(:,:,3) = fphys{1}(:,:,3) - 1/2 * obj.dt .* ( NqHy );

% fphys{1}(:,:,6) = fphys{1}(:,:,6) + obj.dt .* NonhydrostaticPressure;
% fphys{1}(:,:,2) = fphys{1}(:,:,2) - 1/2 * obj.dt .* ( NqHx );
% fphys{1}(:,:,3) = fphys{1}(:,:,3) - 1/2 * obj.dt .* ( NqHy );

end

function [ NqHx , NqHy ] = getNonhydroPartialCorrectionTerm( obj, mesh, WaterDepth, NonhydroPressure )
TempWaterDepth = obj.Vq{1} * WaterDepth;
TempNonhydroPressure  = obj.Vq{1} * NonhydroPressure;

NqHx = ...
    - obj.Dr{1} * ( obj.rxwJ{1} .* ( TempWaterDepth .* TempNonhydroPressure) ) ...
    - obj.Ds{1} * ( obj.sxwJ{1} .* ( TempWaterDepth .* TempNonhydroPressure) ) ...
    - obj.Dt{1} * ( obj.txwJ{1} .* ( TempWaterDepth .* TempNonhydroPressure) );

NqHy =  ...
    - obj.Dr{1} * ( obj.rywJ{1} .* ( TempWaterDepth .* TempNonhydroPressure) ) ...
    - obj.Ds{1} * ( obj.sywJ{1} .* ( TempWaterDepth .* TempNonhydroPressure) ) ...
    - obj.Dt{1} * ( obj.tywJ{1} .* ( TempWaterDepth .* TempNonhydroPressure) );

InnerEdge = mesh.InnerEdge;
[hm, hp] = InnerEdge.matEvaluateSurfValue( num2cell(WaterDepth,[1 2]) );

[Nqm, Nqp] = InnerEdge.matEvaluateSurfValue( num2cell(NonhydroPressure,[1 2]) );

[Ghm, Ghp] = obj.matInterpolateToFaceGaussQuadraturePoint( InnerEdge, obj.IEFVfq{1}, hm, hp);
[GNqm, GNqp] = obj.matInterpolateToFaceGaussQuadraturePoint( InnerEdge, obj.IEFVfq{1}, Nqm, Nqp);

fluxSHX = obj.IEnx{1}.*( Ghm .* GNqm + Ghp .* GNqp )./2;fluxSHY = obj.IEny{1}.*( Ghm .* GNqm + Ghp .* GNqp )./2;

TempIntegralHQX = + ( obj.IELIFT{1} * ( obj.IEwJs{1} .* ( fluxSHX ) ));
NqHx = obj.matAssembleIntoRHS( InnerEdge, TempIntegralHQX, NqHx);
TempIntegralHQY = + ( obj.IELIFT{1} * ( obj.IEwJs{1} .* ( fluxSHY ) ));
NqHy = obj.matAssembleIntoRHS( InnerEdge, TempIntegralHQY, NqHy);

BoundaryEdge = mesh.BoundaryEdge;
[hm, hp] = BoundaryEdge.matEvaluateSurfValue( num2cell(WaterDepth,[1 2]) );

[Nqm, Nqp] = BoundaryEdge.matEvaluateSurfValue( num2cell(NonhydroPressure,[1 2]) );

[Ghm, Ghp] = obj.matInterpolateToFaceGaussQuadraturePoint( BoundaryEdge, obj.BEFVfq{1}, hm, hp);
[GNqm, GNqp] = obj.matInterpolateToFaceGaussQuadraturePoint( BoundaryEdge, obj.BEFVfq{1}, Nqm, Nqp);

fluxSHX = obj.BEnx{1}.*( Ghm .* GNqm + Ghp .* GNqp )./2;fluxSHY = obj.BEny{1}.*( Ghm .* GNqm + Ghp .* GNqp )./2;

TempIntegralX = + ( obj.BELIFT{1} * ( obj.BEwJs{1} .* ( fluxSHX ) ));
TempIntegralY= + ( obj.BELIFT{1} * ( obj.BEwJs{1} .* ( fluxSHY ) ));

NqHx = obj.matAssembleBoundaryAndSourceTermIntoRHS( BoundaryEdge, TempIntegralX, NqHx);
NqHy = obj.matAssembleBoundaryAndSourceTermIntoRHS( BoundaryEdge, TempIntegralY, NqHy);

NqHx = permute( sum( ...
    bsxfun(@times, obj.invM{1}, ...
    permute( permute( NqHx, [1,3,2] ), ...
    [2,1,3] ) ), 2 ), [1,3,2]);

NqHy = permute( sum( ...
    bsxfun(@times, obj.invM{1}, ...
    permute( permute( NqHy, [1,3,2] ), ...
    [2,1,3] ) ), 2 ), [1,3,2]);
end