function [StiffMatrix, NonhydrostaticRHS] = matResembleGlobalMatrix(obj, mesh, StiffMatrix, NonhydrostaticRHS)
%> @brief Function to resemble global stiffMatrix and Nonhydrostatic right hand side
%> @details
%> Function to get the flux term at the wet dry interface
%> @param[in] mesh The mesh object
%> @param[in] FluxTerm The flux term needs to be corrected
%> @param[in] mesh The mesh object
%> @param[out] FluxTerm The flux term corrected

Index = find(mesh.status ~= enumSWERegion.Wet);
data = zeros(mesh.cell.Np, numel(Index));
for i = 1:numel(Index)
    data(:,i) = (Index - 1) * mesh.cell.Np + 1 : Index * mesh.cell.Np;
end
StiffMatrix(data,:) = [];
StiffMatrix(:,data) = [];
NonhydrostaticRHS(data) = [];
end