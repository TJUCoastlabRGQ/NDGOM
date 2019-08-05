function  matReconstructStiffmatrixRelatedMatrix(obj, mesh)
%> @brief Function to assemble the characteristic matrix
%> @details
%> Function to assemble the characteristic matrix with wetting and drying considered
%> @param[in] PhysClass The hydrostatic solver

for i = 1:size(obj.WetDryCell, 2)
    ele = obj.WetDryCell(1 , i );
    [tempPNPX, tempSecondOrderTerm]...
        = obj.matAssembleCharacteristicMatrix( mesh, ele, obj.WetDryCell(2:size(obj.WetDryCell, 1), i)); 
    [obj.TempPNPX(:,(ele-1)* mesh.cell.Np+1:ele * mesh.cell.Np),...
        obj.TempSecondOrderTerm(:,(ele-1) * mesh.cell.Np + 1 : ele * mesh.cell.Np)] = ...
        VectorConvert(tempPNPX, tempSecondOrderTerm);
end

end

function [PNPX, SecondOrderTerm]=VectorConvert(tempPNPX, tempSecondOrderTerm)
PNPX = sparse(tempPNPX);
SecondOrderTerm = sparse(tempSecondOrderTerm);
end