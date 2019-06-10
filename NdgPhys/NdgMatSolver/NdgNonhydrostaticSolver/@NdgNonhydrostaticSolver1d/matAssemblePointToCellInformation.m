function matAssemblePointToCellInformation(obj, K, Np, PNPX, SPNPX, NP )   
tempdata = abs(PNPX) + abs(SPNPX) + NP;
[obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix] = mxAssemblePointToCellInformation(K*Np, tempdata );
end