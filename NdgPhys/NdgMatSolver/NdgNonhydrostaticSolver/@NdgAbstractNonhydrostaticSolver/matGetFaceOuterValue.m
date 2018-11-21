function OuterValue = matGetFaceOuterValue(obj, mesh, OuterValue, InnerValue)

% %> bug exist if N = 3 for quad
OuterValue = mxGetFaceOuterValue(obj.ZeroFluxBoundary, OuterValue, ...
    InnerValue, mesh.cell.TNfp, mesh.cell.Np, mesh.cell.Nfp(1), mesh.cell.Fmask, obj.AdjacentDryCellAndFace);

end