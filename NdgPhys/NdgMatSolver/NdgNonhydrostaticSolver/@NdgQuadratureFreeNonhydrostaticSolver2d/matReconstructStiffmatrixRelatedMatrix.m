function [UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY, UpdatedFNPBX, UpdatedFNPBY]...
    = matReconstructStiffmatrixRelatedMatrix(obj, PhysClass)
%> @brief Function to assemble the characteristic matrix
%> @details
%> Function to assemble the characteristic matrix with wetting and drying considered
%> @param[in] PhysClass The hydrostatic solver
%> @param[out] UpdatedPNPX The PNPX with wetting and drying considered
%> @param[out] UpdatedPNPY The PNPY with wetting and drying considered
%> @param[out] UpdatedSPNPX The SPNPX with wetting and drying considered
%> @param[out] UpdatedSPNPY The SPNPY with wetting and drying considered
%> @param[out] UpdatedFNPBX The UpdatedFNPBX with wetting and drying considered
%> @param[out] UpdatedFNPBY The UpdatedFNPBY with wetting and drying considered
mesh = PhysClass.meshUnion(1);

UpdatedPNPX = obj.TempPNPX; UpdatedPNPY = obj.TempPNPY;
UpdatedSPNPX = obj.TempSPNPX; UpdatedSPNPY = obj.TempSPNPY;

DryToWetPoint = setdiff(obj.TempWetDryPoint, obj.WetDryPoint);

UpdatedPNPX(:,DryToWetPoint) = obj.PNPX(:,DryToWetPoint);
UpdatedPNPY(:,DryToWetPoint) = obj.PNPY(:,DryToWetPoint);
UpdatedSPNPX(:,DryToWetPoint) = obj.SPNPX(:,DryToWetPoint);
UpdatedSPNPY(:,DryToWetPoint) = obj.SPNPY(:,DryToWetPoint);


NewWetDryFace  = setdiff(obj.ZeroFluxBoundary, obj.TempZeroFluxBoundary, 'rows');
TempWetToDryPoint = obj.matGetTempWetToDryPoint(mesh.cell.Np, mesh.cell.Nfp(1), NewWetDryFace, mesh.cell.Fmask);

[obj.NonhydroFmPoint, obj.NonhydroFpPoint] = obj.matAssemblePointRelatedInformation...
    (obj.ZeroFluxBoundary, obj.AdjacentDryCellAndFace, PhysClass.meshUnion(1).InnerEdge.FToE...
    ,PhysClass.meshUnion(1).InnerEdge.FToN1);

if obj.WetNum == mesh.K
    %doing nothing
else
    for i = 1:numel(TempWetToDryPoint)
        num = TempWetToDryPoint(i);
        [tempPNPX, tempPNPY, tempSPNPX, tempSPNPY]...
            = obj.matAssembleCharacteristicMatrix( mesh, num);
        [UpdatedPNPX(:,num), UpdatedPNPY(:,num), UpdatedSPNPX(:,num), UpdatedSPNPY(:,num)]...
            = VectorConvert(tempPNPX, tempPNPY, tempSPNPX, tempSPNPY );
    end
end
obj.TempPNPX = UpdatedPNPX; obj.TempPNPY = UpdatedPNPY;
obj.TempSPNPX = UpdatedSPNPX;  obj.TempSPNPY = UpdatedSPNPY;
obj.TempWetDryPoint = obj.WetDryPoint;
obj.TempZeroFluxBoundary = obj.ZeroFluxBoundary;
end

function [PNPX, PNPY, SPNPX, SPNPY ]=VectorConvert(tempPNPX, tempPNPY, tempSPNPX, tempSPNPY)
PNPX = sparse(tempPNPX);
PNPY = sparse(tempPNPY);
SPNPX = sparse(tempSPNPX);
SPNPY = sparse(tempSPNPY);
end