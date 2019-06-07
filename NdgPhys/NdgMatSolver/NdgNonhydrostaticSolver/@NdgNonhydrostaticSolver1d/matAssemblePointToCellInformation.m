function matAssemblePointToCellInformation(obj, K, Np, PNPX, SPNPX, NP )   
tempdata = abs(PNPX) + abs(PNPY) + abs(SPNPX) + abs(SPNPY) + NP;
[obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix] = mxAssemblePointToCellInformation(K*Np, tempdata );
end