function matAssemblePointToCellInformation(obj, K, Np, PNPX, PNPY, SecondOrderTerm, NP )   
tempdata = abs(PNPX) + abs(PNPY) + abs(SecondOrderTerm) + NP;
[obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix] = mxAssemblePointToCellInformation(K*Np, tempdata );
end