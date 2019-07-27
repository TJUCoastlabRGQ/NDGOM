function matSetInitializeCharacteristicMatrix(obj, physClass, mesh)
%> @brief Function to set the characteristic matrix at the initial stage
%> @details
%> Function to  set the characteristic matrix used when assemble the global stiff matrix
%> @param[in] mesh The mesh object
Num = numel(mesh.x);
K = mesh.K; Np = mesh.cell.Np;
obj.PNPX = spalloc(Num,Num,2*Num*Np); obj.PNPY = spalloc(Num,Num,2*Num*Np);
obj.SecondOrderTerm = spalloc(Num,Num,2*Num*Np); 
obj.NP = spalloc(Num,Num,Num);
for ele =  1:K
    obj.NP( (ele - 1)*Np + 1 : (ele - 1)*Np + 1,(ele - 1)*Np + 1 : (ele - 1)*Np + 1 ) = 1;
    [tempPNPX, tempPNPY, tempSecondOrderTerm]...
        = obj.matAssembleCharacteristicMatrix( mesh, ele, obj.edgeType(:, ele));
    [obj.PNPX(:,(ele-1)*Np+1:ele*Np), obj.PNPY(:,(ele-1)*Np+1:ele*Np),...
        obj.SecondOrderTerm(:,(ele-1)*Np+1:ele*Np)] = ...
        VectorConvert(tempPNPX, tempPNPY, tempSecondOrderTerm);
end

obj.TempPNPX = obj.PNPX; obj.TempPNPY = obj.PNPY;
obj.TempSecondOrderTerm = obj.SecondOrderTerm;
end

function [PNPX, PNPY, SecondOrderTerm]=VectorConvert(tempPNPX, tempPNPY, tempSecondOrderTerm)
PNPX = sparse(tempPNPX);
PNPY = sparse(tempPNPY);
SecondOrderTerm = sparse(tempSecondOrderTerm);
end