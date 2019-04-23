function matAssemblePointToCellInformation(obj, K, Np, PNPX,PNPY, SPNPX, SPNPY, NPBX, NPBY, FNPBX, FNPBY, NP )   
tempdata = abs(PNPX) + abs(PNPY) + abs(SPNPX) + abs(SPNPY) + abs(NPBX) + abs(NPBY) + abs(FNPBX) + abs(FNPBY) + NP;
[obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix] = mxAssemblePointToCellInformation(K*Np, tempdata );
end