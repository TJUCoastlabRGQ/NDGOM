function matSetInitializeCharacteristicMatrix(obj, physClass, mesh)
%> @brief Function to set the characteristic matrix at the initial stage
%> @details
%> Function to  set the characteristic matrix used when assemble the global stiff matrix
%> @param[in] mesh The mesh object
Num = numel(mesh.x);
K = mesh.K; Np = mesh.cell.Np;
obj.PNPX = spalloc(Num,Num,2*Num*Np);
obj.SPNPX = spalloc(Num,Num,2*Num*Np);
obj.NP = spalloc(Num,Num,Num);
for i =  1:K*Np
    [tempPNPX,  tempSPNPX, tempNp]...
        = obj.matAssembleCharacteristicMatrix( mesh, i);
    [obj.PNPX(:,i), obj.SPNPX(:,i), obj.NP(:,i)]...
        = VectorConvert(tempPNPX, tempSPNPX, tempNp);
end

obj.TempPNPX = obj.PNPX;
obj.TempSPNPX = obj.SPNPX;
end

function [PNPX, SPNPX, Np]=VectorConvert(tempPNPX, tempSPNPX, tempNp)
PNPX = sparse(tempPNPX);
SPNPX = sparse(tempSPNPX);
Np = sparse(tempNp);
end