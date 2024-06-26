function  matReconstructStiffmatrixRelatedMatrix(obj, mesh)
%> @brief Function to assemble the characteristic matrix
%> @details
%> Function to assemble the characteristic matrix with wetting and drying considered
%> @param[in] PhysClass The hydrostatic solver

for i = 1:size(obj.WetDryCell, 2)
    ele = obj.WetDryCell(1 , i );
    Np = mesh.cell.Np;
    [tempPNPX, tempPNPY, tempSecondOrderTerm]...
        = obj.matAssembleCharacteristicMatrix( mesh, ele, obj.WetDryCell(2:size(obj.WetDryCell, 1), i));   % abstract in nonhydrostatic solver 2d and 1d, different for Gauss quad and quad free version, placed in the abstract solver;
    [obj.TempPNPX(:,(ele-1)*Np+1:ele*Np), obj.TempPNPY(:,(ele-1)*Np+1:ele*Np),...
        obj.TempSecondOrderTerm(:,(ele-1)*Np+1:ele*Np)] = ...
        VectorConvert(tempPNPX, tempPNPY, tempSecondOrderTerm);
end

end

function [PNPX, PNPY, SecondOrderTerm]=VectorConvert(tempPNPX, tempPNPY, tempSecondOrderTerm)
PNPX = sparse(tempPNPX);
PNPY = sparse(tempPNPY);
SecondOrderTerm = sparse(tempSecondOrderTerm);
end