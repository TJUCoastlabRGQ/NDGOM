function [UpdatedPNPX, UpdatedSPNPX ]= matReconstructStiffmatrixRelatedMatrix(obj, PhysClass)
%> @brief Function to assemble the characteristic matrix
%> @details
%> Function to assemble the characteristic matrix with wetting and drying considered
%> @param[in] PhysClass The hydrostatic solver
%> @param[out] UpdatedPNPX The PNPX with wetting and drying considered
%> @param[out] UpdatedSPNPX The SPNPX with wetting and drying considered
mesh = PhysClass.meshUnion(1);

UpdatedPNPX = obj.TempPNPX;
UpdatedSPNPX = obj.TempSPNPX;

DryToWetPoint = setdiff(obj.TempWetDryPoint, obj.WetDryPoint);

UpdatedPNPX(:,DryToWetPoint) = obj.PNPX(:,DryToWetPoint);
UpdatedSPNPX(:,DryToWetPoint) = obj.SPNPX(:,DryToWetPoint);


NewWetDryFace  = setdiff(obj.ZeroFluxBoundary, obj.TempZeroFluxBoundary, 'rows');
TempWetToDryPoint = obj.matGetTempWetToDryPoint(mesh.cell.Np, mesh.cell.Nfp(1), NewWetDryFace, mesh.cell.Fmask);


if obj.WetNum == mesh.K
    %doing nothing
else
    for i = 1:numel(TempWetToDryPoint)
    [tempPNPX,  tempSPNPX, tempNp]...
        = obj.matAssembleCharacteristicMatrix( mesh, i);
    [obj.PNPX(:,i), obj.SPNPX(:,i), obj.NP(:,i)]...
        = VectorConvert(tempPNPX, tempSPNPX, tempNp);
    end
end
obj.TempPNPX = UpdatedPNPX; 
obj.TempSPNPX = UpdatedSPNPX;
obj.TempWetDryPoint = obj.WetDryPoint;
obj.TempZeroFluxBoundary = obj.ZeroFluxBoundary;
end

function [PNPX, SPNPX , Np]=VectorConvert(tempPNPX, tempSPNPX)
PNPX = sparse(tempPNPX);
SPNPX = sparse(tempSPNPX);
Np = sparse(tempNp);
end