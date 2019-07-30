function matAssemblePointToCellInformation(obj, K, Np )   
tempdata = abs(obj.TempPNPX)  + abs(obj.TempSecondOrderTerm) + obj.NP;
[obj.JcsGlobalStiffMatrix, obj.JrsGlobalStiffMatrix] = mxAssemblePointToCellInformation(K*Np, tempdata );
end