function matSetInitializeCharacteristicMatrix(obj, physClass, mesh)
%> @brief Function to set the characteristic matrix at the initial stage
%> @details
%> Function to  set the characteristic matrix used when assemble the global stiff matrix
%> @param[in] mesh The mesh object
Num = numel(mesh.x);
K = mesh.K; Np = mesh.cell.Np;
obj.PNPX = spalloc(Num,Num,2*Num*Np);obj.PNPY = spalloc(Num,Num,2*Num*Np);
obj.SPNPX = spalloc(Num,Num,2*Num*Np);obj.SPNPY = spalloc(Num,Num,2*Num*Np);
obj.NPBX = spalloc(Num,Num,2*Num*Np);obj.NPBY = spalloc(Num,Num,2*Num*Np);
obj.FNPBX = spalloc(Num,Num,2*Num*Np);obj.FNPBY = spalloc(Num,Num,2*Num*Np);
obj.NP = spalloc(Num,Num,Num);
for i =  1:K*Np
    [tempPNPX, tempPNPY, tempSPNPX, tempSPNPY,tempNPBX, tempNPBY,...
        tempFNPBX, tempFNPBY, tempNp] = obj.matAssembleCharacteristicMatrix( mesh, i);
    [obj.PNPX(:,i), obj.PNPY(:,i), obj.SPNPX(:,i), obj.SPNPY(:,i),...
        obj.NPBX(:,i),obj.NPBY(:,i),obj.FNPBX(:,i), obj.FNPBY(:,i), obj.NP(:,i)] = VectorConvert(tempPNPX, tempPNPY, ...
        tempSPNPX, tempSPNPY,tempNPBX, tempNPBY, tempFNPBX, tempFNPBY, tempNp);
end

obj.TempPNPX = obj.PNPX;
obj.TempPNPY = obj.PNPY;
obj.TempSPNPX = obj.SPNPX;
obj.TempSPNPY = obj.SPNPY;
obj.TempFNPBX = obj.FNPBX;
obj.TempFNPBY = obj.FNPBY;
end

function [PNPX, PNPY, SPNPX, SPNPY, NPBX, NPBY, FNPBX, FNPBY, Np]=VectorConvert(tempPNPX, tempPNPY,...
    tempSPNPX, tempSPNPY,tempNPBX, tempNPBY,tempFNPBX, tempFNPBY, tempNp)
PNPX = sparse(tempPNPX);
PNPY = sparse(tempPNPY);
SPNPX = sparse(tempSPNPX);
SPNPY = sparse(tempSPNPY);
NPBX = sparse(tempNPBX);
NPBY = sparse(tempNPBY);
FNPBX = sparse(tempFNPBX);
FNPBY = sparse(tempFNPBY);
Np = sparse(tempNp);
end