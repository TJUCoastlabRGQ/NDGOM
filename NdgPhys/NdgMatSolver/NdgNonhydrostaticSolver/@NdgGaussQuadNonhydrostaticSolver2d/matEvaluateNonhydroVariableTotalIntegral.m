function [ VolumeIntegralX, VolumeIntegralY ] = matEvaluateNonhydroVariableTotalIntegral( obj, mesh, ftype, VariableX, VariableY, VolumeIntegralX, VolumeIntegralY )
InnerEdge = mesh.InnerEdge;
[fmx, fpx] = InnerEdge.matEvaluateSurfValue( VariableX ); [fmy, fpy] = InnerEdge.matEvaluateSurfValue( VariableY );
[fmx, fpx] = obj.matGetFaceValue(fmx, fpx, ftype); [fmy, fpy] = obj.matGetFaceValue(fmy, fpy, ftype);
[Gfmx, Gfpx] = obj.matInterpolateToFaceGaussQuadraturePoint( InnerEdge, obj.IEFVfq{1}, fmx, fpx);
[Gfmy, Gfpy] = obj.matInterpolateToFaceGaussQuadraturePoint( InnerEdge, obj.IEFVfq{1}, fmy, fpy);

fluxSX = obj.IEnx{1}.*(Gfmx + Gfpx)./2;fluxSY = obj.IEny{1}.*(Gfmy + Gfpy)./2;
TempIntegralX = + ( obj.IELIFT{1} * ( obj.IEwJs{1} .* ( fluxSX ) ));
TempIntegralY = + ( obj.IELIFT{1} * ( obj.IEwJs{1} .* ( fluxSY ) ));

VolumeIntegralX = obj.matAssembleIntoRHS( InnerEdge, TempIntegralX, VolumeIntegralX);
VolumeIntegralY = obj.matAssembleIntoRHS( InnerEdge, TempIntegralY, VolumeIntegralY);



BoundaryEdge = mesh.BoundaryEdge;
[fmy, fpy] = BoundaryEdge.matEvaluateSurfValue( VariableY );    [fmx, fpx] = BoundaryEdge.matEvaluateSurfValue( VariableX );
fpy = obj.matImposeNonhydroRelatedBoundaryCondition(fmy, fpy, ftype, obj.EidBoundaryType);
fpx = obj.matImposeNonhydroRelatedBoundaryCondition(fmx, fpx, ftype, obj.EidBoundaryType);

[Gfmx, Gfpx] = obj.matInterpolateToFaceGaussQuadraturePoint( BoundaryEdge, obj.BEFVfq{1}, fmx, fpx);
[Gfmy, Gfpy] = obj.matInterpolateToFaceGaussQuadraturePoint( BoundaryEdge, obj.BEFVfq{1}, fmy, fpy);

fluxSX = obj.BEnx{1}.*(Gfmx + Gfpx)./2;fluxSY = obj.BEny{1}.*(Gfmy + Gfpy)./2;
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