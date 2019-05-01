function [ UpdatedPNPX, UpdatedPNPY, UpdatedSPNPX, UpdatedSPNPY, UpdatedNP ]...
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
UpdatedNP = obj.NP;
% UpdatedNPBX = obj.NPBX; UpdatedNPBY = obj.NPBY;
% UpdatedFNPBX = obj.TempFNPBX; UpdatedFNPBY = obj.TempFNPBY;

DryToWetPoint = setdiff(obj.TempWetDryPoint, obj.WetDryPoint);

UpdatedPNPX(:,DryToWetPoint) = obj.PNPX(:,DryToWetPoint);
UpdatedPNPY(:,DryToWetPoint) = obj.PNPY(:,DryToWetPoint);
UpdatedSPNPX(:,DryToWetPoint) = obj.SPNPX(:,DryToWetPoint);
UpdatedSPNPY(:,DryToWetPoint) = obj.SPNPY(:,DryToWetPoint);
% UpdatedFNPBX(:,DryToWetPoint) = obj.FNPBX(:,DryToWetPoint);
% UpdatedFNPBY(:,DryToWetPoint) = obj.FNPBY(:,DryToWetPoint);

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
        [ tempNp, tempPNPX, tempPNPY, tempSPNPX, tempSPNPY ] = obj.matAssembleCharacteristicMatrix( mesh, num);
        [UpdatedPNPX(:,num), UpdatedPNPY(:,num), UpdatedSPNPX(:,num), UpdatedSPNPY(:,num), UpdatedNP(:,num)]...
            = VectorConvert(tempPNPX, tempPNPY, tempSPNPX, tempSPNPY,  tempNp);
    end
end
obj.TempPNPX = UpdatedPNPX; obj.TempPNPY = UpdatedPNPY;
obj.TempSPNPX = UpdatedSPNPX;  obj.TempSPNPY = UpdatedSPNPY;
% obj.TempFNPBX = UpdatedFNPBX; obj.TempFNPBY = UpdatedFNPBY;
obj.NP = UpdatedNP;
obj.TempWetDryPoint = obj.WetDryPoint;
obj.TempZeroFluxBoundary = obj.ZeroFluxBoundary;
end

function [PNPX, PNPY, SPNPX, SPNPY, Np]=VectorConvert(tempPNPX, tempPNPY,...
    tempSPNPX, tempSPNPY, tempNp)
PNPX = sparse(tempPNPX);
PNPY = sparse(tempPNPY);
SPNPX = sparse(tempSPNPX);
SPNPY = sparse(tempSPNPY);
% NPBX = sparse(tempNPBX);
% NPBY = sparse(tempNPBY);
% FNPBX = sparse(tempFNPBX);
% FNPBY = sparse(tempFNPBY);
Np = sparse(tempNp);
end