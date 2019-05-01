function [ VolumeIntegralX, VolumeIntegralY ] = matEvaluateConservativeVariableTotalIntegral( obj, physClass, mesh, ftype, fphys, VolumeIntegralX, VolumeIntegralY, index )
% [ Gfm, Gfp ] = obj.matInterpolateToFaceGaussQuadraturePoint( edge, FVfq{1}, fm, fp);
% FluxX =
% mesh = physClass.meshUnion(1);
InnerEdge = mesh.InnerEdge;
[fm, fp] = InnerEdge.matEvaluateSurfValue( fphys );
[fm, fp] = obj.matGetFaceValue(fm(:,:,index), fp(:,:,index), ftype);
[Gfm, Gfp] = obj.matInterpolateToFaceGaussQuadraturePoint( InnerEdge, obj.IEFVfq{1}, fm, fp);
fluxSX = obj.IEnx{1}.*(Gfm + Gfp)./2;fluxSY = obj.IEny{1}.*(Gfm + Gfp)./2;
TempIntegralX = + ( obj.IELIFT{1} * ( obj.IEwJs{1} .* ( fluxSX ) ));
VolumeIntegralX = obj.matAssembleIntoRHS( InnerEdge, TempIntegralX, VolumeIntegralX);
TempIntegralY = + ( obj.IELIFT{1} * ( obj.IEwJs{1} .* ( fluxSY ) ));
VolumeIntegralY = obj.matAssembleIntoRHS( InnerEdge, TempIntegralY, VolumeIntegralY);

BoundaryEdge = mesh.BoundaryEdge;
[fm, fp] = BoundaryEdge.matEvaluateSurfValue( fphys ); 
[fm, fp] = physClass.matImposeBoundaryCondition( BoundaryEdge, BoundaryEdge.nx, BoundaryEdge.ny, fm, fp, physClass.fext{1} );
[Gfm, Gfp] = obj.matInterpolateToFaceGaussQuadraturePoint( BoundaryEdge, obj.BEFVfq{1}, fm(:,:,index), fp(:,:,index));
fluxSX = obj.BEnx{1}.*(Gfm + Gfp)./2;fluxSY = obj.BEny{1}.*(Gfm + Gfp)./2;
TempIntegralX = + ( obj.BELIFT{1} * ( obj.BEwJs{1} .* ( fluxSX ) ));
TempIntegralY= + ( obj.BELIFT{1} * ( obj.BEwJs{1} .* ( fluxSY ) ));

VolumeIntegralX = obj.matAssembleBoundaryAndSourceTermIntoRHS( BoundaryEdge, TempIntegralX, VolumeIntegralX);
VolumeIntegralY = obj.matAssembleBoundaryAndSourceTermIntoRHS( BoundaryEdge, TempIntegralY, VolumeIntegralY);

VolumeIntegralX = permute( sum( ...
    bsxfun(@times, obj.invM{1}, ...
    permute( permute( VolumeIntegralX, [1,3,2] ), ...
    [2,1,3] ) ), 2 ), [1,3,2]);

VolumeIntegralY = permute( sum( ...
    bsxfun(@times, obj.invM{1}, ...
    permute( permute( VolumeIntegralY, [1,3,2] ), ...
    [2,1,3] ) ), 2 ), [1,3,2]);
end