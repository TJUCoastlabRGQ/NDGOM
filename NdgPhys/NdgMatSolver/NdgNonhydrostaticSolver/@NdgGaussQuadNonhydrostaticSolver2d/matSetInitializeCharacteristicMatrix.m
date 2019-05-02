function matSetInitializeCharacteristicMatrix(obj, physClass, mesh)
%> @brief Function to set the characteristic matrix at the initial stage
%> @details
%> Function to  set the characteristic matrix used when assemble the global stiff matrix
%> @param[in] mesh The mesh object
Num = numel(obj.wJ{1});
K = mesh.K; Np = mesh.cell.Np;
obj.PNPX = spalloc(Num,K*Np,2*Num*Np);obj.PNPY = spalloc(Num,K*Np,2*Num*Np);
obj.SPNPX = spalloc(Num,K*Np,2*Num*Np);obj.SPNPY = spalloc(Num,K*Np,2*Num*Np);
% obj.NPBX = spalloc(Num,Num,2*Num*Np);obj.NPBY = spalloc(Num,Num,2*Num*Np);
% obj.FNPBX = spalloc(Num,Num,2*Num*Np);obj.FNPBY = spalloc(Num,Num,2*Num*Np);
obj.NP = spalloc(Num,K*Np,Num);
for i =  1:K*Np
    [ tempNp, tempPNPX, tempPNPY, tempSPNPX, tempSPNPY ] = obj.matAssembleCharacteristicMatrix( mesh, i);
    [ obj.PNPX(:,i), obj.PNPY(:,i), obj.SPNPX(:,i), obj.SPNPY(:,i),...
         obj.NP(:,i) ] = VectorConvert(tempPNPX, tempPNPY, tempSPNPX, tempSPNPY, tempNp);
end

obj.TempPNPX = obj.PNPX;
obj.TempPNPY = obj.PNPY;
obj.TempSPNPX = obj.SPNPX;
obj.TempSPNPY = obj.SPNPY;
% obj.TempFNPBX = obj.FNPBX;
% obj.TempFNPBY = obj.FNPBY;
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