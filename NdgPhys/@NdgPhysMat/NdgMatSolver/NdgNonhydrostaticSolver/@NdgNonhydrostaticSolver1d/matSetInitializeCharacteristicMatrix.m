function matSetInitializeCharacteristicMatrix(obj, mesh)
%> @brief Function to set the characteristic matrix at the initial stage
%> @details
%> Function to  set the characteristic matrix used when assemble the global stiff matrix
%> @param[in] mesh The mesh object
Num = numel(mesh.x);
K = mesh.K; Np = mesh.cell.Np;
obj.PNPX = spalloc(Num,Num,2*Num*Np);
obj.SecondOrderTerm = spalloc(Num,Num,2*Num*Np); 
obj.NP = speye(K*Np, K*Np);
for ele =  1:K
    [tempPNPX, tempSecondOrderTerm]...
        = obj.matAssembleCharacteristicMatrix( mesh, ele, obj.edgeType(:, ele));   % abstract in nonhydrostatic solver 2d and 1d, different for Gauss quad and quad free version, placed in the abstract solver;
    [obj.PNPX(:,(ele-1)*Np+1:ele*Np), ...
        obj.SecondOrderTerm(:,(ele-1)*Np+1:ele*Np)] = ...
        VectorConvert(tempPNPX, tempSecondOrderTerm);
end

obj.TempPNPX = obj.PNPX; 
obj.TempSecondOrderTerm = obj.SecondOrderTerm;
end

function [PNPX, SecondOrderTerm]=VectorConvert(tempPNPX, tempSecondOrderTerm)
PNPX = sparse(tempPNPX);
SecondOrderTerm = sparse(tempSecondOrderTerm);
end