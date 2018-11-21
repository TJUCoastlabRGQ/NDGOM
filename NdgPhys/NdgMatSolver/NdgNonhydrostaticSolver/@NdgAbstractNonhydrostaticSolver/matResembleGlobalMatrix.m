function [StiffMatrix, NonhydrostaticRHS] = matResembleGlobalMatrix(obj, mesh, StiffMatrix, NonhydrostaticRHS)
%> @brief Function to resemble global stiffMatrix and Nonhydrostatic right hand side
%> @details
%> Function to get the flux term at the wet dry interface
%> @param[in] mesh The mesh object
%> @param[in] FluxTerm The flux term needs to be corrected
%> @param[in] mesh The mesh object
%> @param[out] FluxTerm The flux term corrected
if obj.WetNum == mesh.K
    %doing nothing
else
    cellIndex = 1:mesh.K;
    DryCell = setdiff(cellIndex, obj.WetCellIndex);
    data = zeros(mesh.cell.Np, numel(DryCell));
    for i = 1:numel(DryCell)
        data(:,i) = ((DryCell(i) - 1)* mesh.cell.Np+1 :DryCell(i)* mesh.cell.Np)';
    end
    StiffMatrix(data,:) = [];
    StiffMatrix(:,data) = [];
    NonhydrostaticRHS(data) = [];
end
end