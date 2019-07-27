function  matReconstructStiffmatrixRelatedMatrix(obj)
%> @brief Function to assemble the characteristic matrix
%> @details
%> Function to assemble the characteristic matrix with wetting and drying considered
%> @param[in] PhysClass The hydrostatic solver

for i = 1:size(obj.WetDryCell, 2)
    ele = obj.WetDryCell(1 , i );
    [tempPNPX, tempPNPY, tempSecondOrderTerm]...
        = obj.matAssembleCharacteristicMatrix( mesh, ele, obj.WetDryCell(2:size(obj.WetDryCell, 1), i));
    [obj.TempPNPX(:,(ele-1)*Np+1:ele*Np), obj.TempPNPY(:,(ele-1)*Np+1:ele*Np),...
        obj.TempSecondOrderTerm(:,(ele-1)*Np+1:ele*Np)] = ...
        VectorConvert(tempPNPX, tempPNPY, tempSecondOrderTerm);
end

end

function [PNPX, PNPY, SPNPX, SPNPY ]=VectorConvert(tempPNPX, tempPNPY, tempSPNPX, tempSPNPY)
PNPX = sparse(tempPNPX);
PNPY = sparse(tempPNPY);
SPNPX = sparse(tempSPNPX);
SPNPY = sparse(tempSPNPY);
end