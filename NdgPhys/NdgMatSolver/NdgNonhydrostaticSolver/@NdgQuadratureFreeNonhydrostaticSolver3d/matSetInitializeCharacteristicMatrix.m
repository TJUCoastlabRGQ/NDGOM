function matSetInitializeCharacteristicMatrix( obj, meshUnion )

    obj.SPNPX = mxSecondOrderDerivAboutNohydroPressInHorizon( meshUnion.cell.Np,...
        meshUnion.K, meshUnion.Nz, meshUnion.mesh2d.K, meshUnion.cell.Nface,...
        meshUnion.EToE, meshUnion.InnerEdge.FToE, meshUnion.InnerEdge.FToF,...
        meshUnion.LAV, meshUnion.J, meshUnion.InnerEdge.Js, meshUnion.cell.M,...
        meshUnion.InnerEdge.M, meshUnion.BoundaryEdge.Ne , meshUnion.InnerEdge.Ne, ...
        meshUnion.InnerEdge.LAV, meshUnion.cell.Fmask, meshUnion.cell.Nfp, meshUnion.cell.N,...
        meshUnion.InnerEdge.nx, meshUnion.rx, meshUnion.sx, meshUnion.cell.Dr, meshUnion.cell.Ds);
    
    obj.SPNPY = mxSecondOrderDerivAboutNohydroPressInHorizon( meshUnion.cell.Np,...
        meshUnion.K, meshUnion.Nz, meshUnion.mesh2d.K, meshUnion.cell.Nface,...
        meshUnion.EToE, meshUnion.InnerEdge.FToE, meshUnion.InnerEdge.FToF,...
        meshUnion.LAV, meshUnion.J, meshUnion.InnerEdge.Js, meshUnion.cell.M,...
        meshUnion.InnerEdge.M, meshUnion.BoundaryEdge.Ne , meshUnion.InnerEdge.Ne, ...
        meshUnion.InnerEdge.LAV, meshUnion.cell.Fmask, meshUnion.cell.Nfp, meshUnion.cell.N,...
        meshUnion.InnerEdge.ny, meshUnion.ry, meshUnion.sy, meshUnion.cell.Dr, meshUnion.cell.Ds);    
    
    obj.SPNPS = mxSecondOrderDerivAboutNohydroPressInVert( meshUnion.cell.Np,...
        meshUnion.K, meshUnion.Nz, meshUnion.mesh2d.K, meshUnion.cell.Nface,...
        meshUnion.EToE, meshUnion.tz, meshUnion.cell.Dt, meshUnion.J, meshUnion.mesh2d.J,...
        meshUnion.mesh2d.cell.M, meshUnion.cell.M, meshUnion.cell.Fmask, meshUnion.cell.Nfp,...
        meshUnion.cell.N);

    obj.PNPS = mxFirstOrderDerivAboutNohydroPressInVert( meshUnion.cell.Np,...
        meshUnion.K, meshUnion.Nz, meshUnion.mesh2d.K, meshUnion.cell.Nface,...
        meshUnion.EToE, meshUnion.tz, meshUnion.cell.Dt, meshUnion.J, meshUnion.mesh2d.J,...
        meshUnion.mesh2d.cell.M, meshUnion.cell.M, meshUnion.cell.Fmask, meshUnion.cell.Nfp);
    
    obj.MSPNPX = mxMixedSecondOrderDerivAboutNonhydroStaticPressureUpdated(meshUnion.cell.Np,...
        meshUnion.K, meshUnion.Nz, meshUnion.mesh2d.K, meshUnion.cell.Nface,...
        meshUnion.EToE, meshUnion.InnerEdge.FToE, meshUnion.InnerEdge.FToF,...
        meshUnion.LAV, meshUnion.J, meshUnion.mesh2d.J, meshUnion.InnerEdge.Js,...
        meshUnion.cell.M, meshUnion.InnerEdge.M, meshUnion.mesh2d.cell.M, ...
        meshUnion.BoundaryEdge.Ne , meshUnion.InnerEdge.Ne, meshUnion.InnerEdge.LAV, ...
        meshUnion.cell.Fmask, meshUnion.cell.Nfp, meshUnion.cell.N, meshUnion.InnerEdge.nx,...
        meshUnion.rx, meshUnion.sx, meshUnion.tz, meshUnion.cell.Dr, meshUnion.cell.Ds,...
        meshUnion.cell.Dt, meshUnion.BoundaryEdge.FToF, meshUnion.BoundaryEdge.FToE, ...
        meshUnion.BoundaryEdge.nx, meshUnion.BoundaryEdge.Js, meshUnion.BottomEdge.LAV, ...
        meshUnion.BottomEdge.Ne, meshUnion.BottomEdge.FToE, meshUnion.BottomEdge.FToF,...
        meshUnion.SurfaceBoundaryEdge.LAV, meshUnion.SurfaceBoundaryEdge.Ne, meshUnion.SurfaceBoundaryEdge.FToE);
    
    obj.MSPNPY = mxMixedSecondOrderDerivAboutNonhydroStaticPressureUpdated(meshUnion.cell.Np,...
        meshUnion.K, meshUnion.Nz, meshUnion.mesh2d.K, meshUnion.cell.Nface,...
        meshUnion.EToE, meshUnion.InnerEdge.FToE, meshUnion.InnerEdge.FToF,...
        meshUnion.LAV, meshUnion.J, meshUnion.mesh2d.J, meshUnion.InnerEdge.Js,...
        meshUnion.cell.M, meshUnion.InnerEdge.M, meshUnion.mesh2d.cell.M, ...
        meshUnion.BoundaryEdge.Ne , meshUnion.InnerEdge.Ne, meshUnion.InnerEdge.LAV, ...
        meshUnion.cell.Fmask, meshUnion.cell.Nfp, meshUnion.cell.N, meshUnion.InnerEdge.ny,...
        meshUnion.ry, meshUnion.sy, meshUnion.tz, meshUnion.cell.Dr, meshUnion.cell.Ds,...
        meshUnion.cell.Dt, meshUnion.BoundaryEdge.FToF, meshUnion.BoundaryEdge.FToE, ...
        meshUnion.BoundaryEdge.ny, meshUnion.BoundaryEdge.Js, meshUnion.BottomEdge.LAV, ...
        meshUnion.BottomEdge.Ne, meshUnion.BottomEdge.FToE, meshUnion.BottomEdge.FToF,...
        meshUnion.SurfaceBoundaryEdge.LAV, meshUnion.SurfaceBoundaryEdge.Ne, meshUnion.SurfaceBoundaryEdge.FToE);   
end